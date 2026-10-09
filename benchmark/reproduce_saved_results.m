function result = reproduce_saved_results(inputDir,outputDir,options)
%REPRODUCE_SAVED_RESULTS Analyze a completed capture without solver calls.
%   RESULT = REPRODUCE_SAVED_RESULTS(INPUTDIR,OUTPUTDIR) physically copies the
%   saved raw inputs into fresh OUTPUTDIR/captured_replay, verifies their
%   hashes, replays the captured summary, and publishes selected CSV tables.
%   Original files remain unchanged. No solver, worker pool or timer is run
%   to measure tension distribution performance. Analysis elapsed time is bookkeeping.
%
%   OPTIONS.interior_table=false: optional independent-readback-based interior
%   table, requiring INPUTDIR/independent_readback_audit.mat from the original
%   physical audit. OPTIONS.norm_statistics=false: optional final DM-relative
%   L2 statistics for the original full run, requiring OPTIONS.inventory_file
%   and original masks, T_quality, completion receipt and source snapshots.
%   Norm analysis is deliberately unsupported on the small smoke profile.
%
%   Saved-data replay uses base MATLAB R2024b and Java. Neither Optimization
%   nor Parallel Computing Toolbox is required. Full replay needs about
%   12 GB free disk for physical raw copies and several GB of working RAM.
%   Generated captured_replay/T_L2_matched.csv is the old QP-relative analysis,
%   NOT the current DM-relative L2 result. Only the optional norm_statistics
%   path produces the latter; see README.txt for their distinct denominators.
%
%   Example:
%     addpath('benchmark');
%     reproduce_saved_results(fullfile(tempdir,'review_smoke_new'), ...
%         fullfile(tempdir,'review_smoke_analysis_new'));

assert(nargin>=2,'Provide input and fresh output directories.');
if nargin<3,options=struct();end
defaults=struct('interior_table',false, ...
    'norm_statistics',false,'inventory_file','');
assert(isstruct(options)&&isscalar(options),'Options must be a scalar struct.');
unknown=setdiff(fieldnames(options),fieldnames(defaults));
assert(isempty(unknown),'Unknown options: %s',strjoin(unknown,', '));
names=fieldnames(defaults);
for k=1:numel(names),if ~isfield(options,names{k}),options.(names{k})=defaults.(names{k});end,end
assert(usejava('jvm'),'Java is required for SHA-256.');
inputDir=canonical(inputDir);outputDir=canonical(outputDir);
assert(isfolder(inputDir),'Input directory is missing.');
assert(~isfolder(outputDir)&&~isfile(outputDir),'Choose a fresh output directory.');
assert(~within(inputDir,outputDir)&&~within(outputDir,inputDir), ...
    'Input/output directory trees must be disjoint.');
here=fileparts(mfilename('fullpath'));captured=fullfile(here,'captured');
loaded=load(fullfile(inputDir,'manifest.mat'),'manifest');manifest=loaded.manifest;
assert(any(strcmp(manifest.status,{'raw_complete','complete'})),'Input capture is incomplete.');
if options.norm_statistics
    assert(manifest.ninputs==30965760,'Norm reproduction requires the original full saved run.');
    assert(isfile(options.inventory_file),'Norm reproduction needs the original inventory file.');
    options.inventory_file=canonical(options.inventory_file);
end
files=[{'manifest.mat','geometry.mat','wec.mat'},strcat(manifest.files,'.mat')];
if options.interior_table,files=[files,{'independent_readback_audit.mat'}];end
for k=1:numel(files),assert(isfile(fullfile(inputDir,files{k})),'Required input missing: %s',files{k});end
mkdir(outputDir);replay=fullfile(outputDir,'captured_replay');mkdir(replay);
records=repmat(struct('file','','sha256_before','','sha256_copy','','sha256_after',''),numel(files),1);
for k=1:numel(files)
    src=fullfile(inputDir,files{k});dst=fullfile(replay,files{k});
    records(k).file=files{k};records(k).sha256_before=sha256(src);
    copyfile(src,dst);records(k).sha256_copy=sha256(dst);
    assert(strcmp(records(k).sha256_before,records(k).sha256_copy),'Copy differs: %s',files{k});
    fprintf('Copied verified saved input: %s\n',files{k});
end
guardDir=fullfile(outputDir,'solver_guards');mkdir(guardDir);
for k=1:numel(manifest.solver_functions)
    fn=manifest.solver_functions{k};fid=fopen(fullfile(guardDir,[fn '.m']),'w');assert(fid>=0);
    fprintf(fid,'function varargout=%s(varargin)\nerror(''saved_replay:SolverForbidden'',''Solver calls are forbidden during saved-data reproduction.'');\nend\n',fn);
    fclose(fid);
end
oldPath=path;oldFolder=pwd;cleanup=onCleanup(@()restore_environment(oldPath,oldFolder)); %#ok<NASGU>
addpath(here,'-begin');addpath(captured,'-begin');addpath(guardDir,'-begin');cd(guardDir);
clear summarize_results plot_results analyze_interior_infeasible
clear(manifest.solver_functions{:});rehash;
for k=1:numel(manifest.solver_functions)
    fn=manifest.solver_functions{k};assert(strcmpi(which(fn),fullfile(guardDir,[fn '.m'])));
end
summary=summarize_results(replay,false);
publish={'T_domain.csv','T_timing.csv','T_quality.csv','T_postwec_components.csv', ...
    'T_iterative_domains.csv','T_wec_behavior.csv','T_pipeline_parity.csv', ...
    'T_force_approximation.csv','T_complete_paths_matched.csv'};
if options.interior_table
    analyze_interior_infeasible(replay);
    publish=[publish,{'T_complete_paths_matched_interiors.csv'}];
end
for k=1:numel(publish),copyfile(fullfile(replay,publish{k}),fullfile(outputDir,publish{k}));end
if options.norm_statistics
    compute_norm_differences_saved(inputDir,fullfile(outputDir,'norm_analysis'),options.inventory_file);
end
for k=1:numel(files)
    records(k).sha256_after=sha256(fullfile(inputDir,files{k}));
    assert(strcmp(records(k).sha256_before,records(k).sha256_after),'Original input changed: %s',files{k});
end
sources={'reproduce_saved_results.m','captured/summarize_results.m', ...
    'captured/analyze_interior_infeasible.m'};
sourceHashes=repmat(struct('file','','sha256',''),numel(sources),1);
for k=1:numel(sources),sourceHashes(k)=struct('file',sources{k},'sha256',sha256(fullfile(here,sources{k})));end
result=struct('status','passed','ninputs',manifest.ninputs,'input_files_unchanged',true, ...
    'new_solver_calls',0,'new_timing_observations',0,'published_tables',{publish}, ...
    'options',options,'input_hashes',records,'analysis_sources',sourceHashes, ...
    'warning','captured_replay/T_L2_matched.csv is historical QP-relative output; use norm_analysis for current DM-relative norms');
fid=fopen(fullfile(outputDir,'saved_reproduction_receipt.json'),'w');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(result,'PrettyPrint',true));fclose(fid);
fid=fopen(fullfile(outputDir,'TABLES_README.txt'),'w');assert(fid>=0);
fprintf(fid,'Selected timing, quality, domain and matched complete-path tables were replayed from saved raw observations. All attempted failures and timing tails are retained.\n');
fprintf(fid,'The captured_replay subfolder retains exact historical helper outputs, including a QP-relative L2 table and descriptive report. They are not the final DM-relative L2 analysis or paper figures.\n');
fprintf(fid,'The current DM-relative L2 statistics are generated only by options.norm_statistics=true into norm_analysis, requiring the original full-run input and inventory.\n');
fclose(fid);
fprintf('Saved-data reproduction passed: %d inputs; no solver calls.\n',summary.manifest.ninputs);
end

function value=canonical(value)
value=char(value);
if ~java.io.File(value).isAbsolute(),value=fullfile(pwd,value);end
value=char(java.io.File(value).getCanonicalPath());
end
function yes=within(parent,child)
yes=strcmpi(parent,child)||startsWith(lower(child),[lower(parent) filesep]);
end
function restore_environment(oldPath,oldFolder)
cd(oldFolder);path(oldPath);
end
function value=sha256(filename)
fid=fopen(filename,'rb');assert(fid>=0,'Cannot read %s',filename);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),data=fread(fid,8*1024*1024,'*uint8');md.update(data);end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
