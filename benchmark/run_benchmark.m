function result = run_benchmark(outputDir,options)
%RUN_BENCHMARK Capture a new single-pass benchmark in a fresh directory.
%   RESULT = RUN_BENCHMARK(OUTPUTDIR) runs the 192-input smoke profile using
%   two process workers, without affinity changes or subsequent analysis.
%   RESULT = RUN_BENCHMARK(OUTPUTDIR,OPTIONS) accepts profile ('smoke',
%   'representative', 'full'), workers, pin_pcores, chunk_size and seed.
%   Full explicitly selects 30,965,760 inputs and defaults to eight workers
%   and the captured hardware-conditional affinity policy. Numerical solver
%   settings and measured call boundaries are preserved in the captured code.
%
%   Requires MATLAB R2024b, Parallel Computing Toolbox, Optimization Toolbox,
%   and Java. A newly created worker pool is closed on return; an existing
%   process pool must match OPTIONS.workers and is left open. There is one
%   wall-clock observation per attempted method/input, including failures.
%   Run REPRODUCE_SAVED_RESULTS separately to analyze a completed capture.
%
%   Example, from the extracted package root:
%     addpath('benchmark');
%     run_benchmark(fullfile(tempdir,'review_smoke_new'));

if nargin<2,options=struct();end
assert(nargin>=1&&~isempty(outputDir),'Provide a fresh output directory.');
assert(isstruct(options)&&isscalar(options),'Options must be a scalar struct.');
valid={'profile','workers','pin_pcores','chunk_size','seed'};
unknown=setdiff(fieldnames(options),valid);
assert(isempty(unknown),'Unknown options: %s',strjoin(unknown,', '));
if ~isfield(options,'profile'),options.profile='smoke';end
options.profile=char(options.profile);
assert(ismember(options.profile,{'smoke','representative','full'}),'Unknown profile.');
if ~isfield(options,'workers')
    if strcmp(options.profile,'smoke'),options.workers=2;else,options.workers=8;end
end
if ~isfield(options,'pin_pcores'),options.pin_pcores=~strcmp(options.profile,'smoke');end
if ~isfield(options,'chunk_size'),options.chunk_size=100000;end
if ~isfield(options,'seed'),options.seed=20261006;end
assert(license('test','Distrib_Computing_Toolbox')&&exist('parpool','file')~=0, ...
    'The measured parallel harness requires Parallel Computing Toolbox.');
assert(license('test','Optimization_Toolbox')&&exist('quadprog','file')~=0, ...
    'The bounded-QP comparator requires Optimization Toolbox.');
assert(usejava('jvm'),'Java is required for source hashes.');
here=fileparts(mfilename('fullpath'));sourceDir=fullfile(fileparts(here),'src');
outputDir=char(outputDir);
if ~java.io.File(outputDir).isAbsolute(),outputDir=fullfile(pwd,outputDir);end
outputDir=char(java.io.File(outputDir).getCanonicalPath());
assert(~isfolder(outputDir)&&~isfile(outputDir),'Choose a new output directory.');
functions={'wt4_2024_minmax_v4','wt_pott_v2','wt_qp_bounded_v11', ...
    'wt_gouttefarde_v2','wt_acta','wt_xacta','WEC_v5'};
for k=1:numel(functions)
    assert(isfile(fullfile(sourceDir,[functions{k} '.m'])),'Missing package source: %s',functions{k});
end
mkdir(outputDir);runtimeDir=fullfile(outputDir,'runtime_source');mkdir(runtimeDir);
copyfile(fullfile(here,'captured','Compare_Review_v11.m'),runtimeDir);
for k=1:numel(functions)
    copyfile(fullfile(sourceDir,[functions{k} '.m']),runtimeDir);
end
oldPath=path;oldFolder=pwd;createdPool=isempty(gcp('nocreate'));
cleanup=onCleanup(@()restore_environment(oldPath,oldFolder,createdPool)); %#ok<NASGU>
addpath(runtimeDir,'-begin');cd(runtimeDir);
clear Compare_Review_v11
clear(functions{:});rehash;
assert(strcmpi(which('Compare_Review_v11'),fullfile(runtimeDir,'Compare_Review_v11.m')));
config=options;config.outdir=outputDir;config.analyze=false;
fprintf('New %s capture; results: %s\n',config.profile,outputDir);
fprintf('This run will not select or replace recorded timing observations.\n');
result=Compare_Review_v11(config);
assert(strcmp(result.status,'raw_complete'),'Capture did not complete.');
end

function restore_environment(oldPath,oldFolder,createdPool)
cd(oldFolder);path(oldPath);
if createdPool
    pool=gcp('nocreate');
    if ~isempty(pool),delete(pool);end
end
end
