function result = compute_norm_differences_saved(inputDir,outputDir,inventoryFile)
%COMPUTE_NORM_DIFFERENCES_SAVED Actual-tension L2 differences from saved T.
%   RESULT = COMPUTE_NORM_DIFFERENCES_SAVED(INPUT,OUTPUT,INVENTORY)
%   computes ||T_method||_2 - ||T_DM||_2 in N for ICFM, bounded QP and VTDA-L2.
%   No solver is called, no timing is measured/replaced, and INPUT is read-only.
%   INPUT is the completed review_v11_full directory. INVENTORY defaults to
%   ../../FULL_DATA_INVENTORY.csv relative to INPUT. OUTPUT may already contain
%   this script, but it must not contain earlier norm-difference outputs.
%
%   Start mask: frame interior, checked WEC code 1, actual unchanged-force
%   flag (<=1e-12 N), exactly 24,433,774 rows. Both members of every retained
%   pair must have code 1, four finite actual tensions, bound violation
%   <=1e-9 N and ||W*T-F_raw||_2<=1e-8 N. QP additionally needs its positive
%   native exitflag. DM/ICFM/VTDA expose no convergence flag in this capture.
%   Actual W and raw force are independently reconstructed from the manifest.
%
%   Quantiles: explicit Hyndman-Fan type 5, midpoint positions (i-0.5)/N,
%   linearly interpolated with endpoint clamping. Zero and negative objective
%   differences remain included. Timings are not read or filtered.
%
%   Outputs include statistics CSV, overlapping and mutually exclusive
%   exclusion counts, compressed full-row masks/reason bits, extrema with
%   actual tensions, readback checks, before/after hashes and a JSON receipt.

assert(nargin>=1&&~isempty(inputDir),'Specify the saved full-run directory.');
here=fileparts(mfilename('fullpath'));
if nargin<2||isempty(outputDir),outputDir=here;end
if nargin<3||isempty(inventoryFile)
    inventoryFile=fullfile(fileparts(fileparts(inputDir)),'FULL_DATA_INVENTORY.csv');
end
inputDir=canonical(inputDir);outputDir=canonical(outputDir);
assert(isfolder(inputDir),'Input is missing.');
assert(~same_or_child(inputDir,outputDir)&&~same_or_child(outputDir,inputDir), ...
    'Input/output directory trees must be disjoint.');
if ~isfolder(outputDir),mkdir(outputDir);end
assert(~isfile(fullfile(outputDir,'norm_difference_receipt.json'))&& ...
    ~isfile(fullfile(outputDir,'norm_difference_statistics.csv')), ...
    'Previous outputs exist; choose a fresh output directory.');
started_at=datestr(now,30);
receipt=struct('status','running','started_at',started_at,'input_directory',inputDir, ...
    'output_directory',outputDir,'new_solver_calls',0,'new_timing_observations',0);
old_dir=pwd;old_path=path;cleanup=onCleanup(@()restore_environment(old_dir,old_path)); %#ok<NASGU>
try
    inventory=readtable(inventoryFile,'TextType','string','VariableNamingRule','preserve');
    receipt.inventory_file=canonical(inventoryFile);receipt.inventory_sha256=hash_file(inventoryFile);
    names={'manifest.mat','geometry.mat','wec.mat','masks.mat','T_quality.csv', ...
        'dm_postwec.mat','icfm_postwec.mat','qp_postwec.mat','vtda_postwec.mat', ...
        'completion_receipt.json','source_snapshot/Compare_Review_v11.m', ...
        'source_snapshot/wt4_2024_minmax_v4.m','source_snapshot/wt_pott_v2.m', ...
        'source_snapshot/wt_qp_bounded_v11.m','source_snapshot/wt_gouttefarde_v2.m'};
    fprintf('Verifying saved-run input hashes.\n');
    before=check_files(inputDir,inventory,names);
    writetable(before,fullfile(outputDir,'norm_input_hashes_before.csv'));
    assert(all(before.passed),'An original input differs from its inventory.');
    s=load(fullfile(inputDir,'manifest.mat'),'manifest','WDM');manifest=s.manifest;WDM=s.WDM;
    assert(ismember(manifest.status,{'raw_complete','complete'})&&manifest.ninputs==30965760);
    assert(manifest.config.force_tolerance==1e-8&&manifest.config.bound_tolerance==1e-9);
    assert(manifest.config.unchanged_tolerance==1e-12);
    n=manifest.ninputs;np=manifest.nposes;nf=manifest.nforces;
    assert(np==4096&&nf==7560&&n==np*nf);
    assert(all(WDM.param.lim_inf(:)==3.8)&&all(WDM.param.lim_sup(:)==25));
    method_names={'ICFM','Bounded QP','VTDA-L2'};
    file_names={'dm_postwec','icfm_postwec','qp_postwec','vtda_postwec'};
    native_applicable=[false false true false];

    sourceFile=[mfilename('fullpath') '.m'];sourceHash=hash_file(sourceFile);
    sourceDir=fullfile(outputDir,'norm_difference_source');mkdir(sourceDir);
    copyfile(sourceFile,fullfile(sourceDir,'compute_norm_differences_saved.m'));
    guardDir=fullfile(outputDir,'norm_solver_guards');mkdir(guardDir);
    for k=1:numel(manifest.solver_functions)
        fn=manifest.solver_functions{k};fid=fopen(fullfile(guardDir,[fn '.m']),'w');assert(fid~=-1);
        fprintf(fid,'function varargout=%s(varargin)\nerror(''saved_norm:SolverForbidden'',''Solver calls are forbidden in this saved-data postprocessor.'');\nend\n',fn);
        fclose(fid);
    end
    addpath(guardDir,'-begin');cd(guardDir);clear(manifest.solver_functions{:});rehash;
    for k=1:numel(manifest.solver_functions)
        fn=manifest.solver_functions{k};assert(strcmpi(which(fn),fullfile(guardDir,[fn '.m'])));
    end

    w=load(fullfile(inputDir,'wec.mat'),'code','unchanged','completed_chunks');
    assert(all(w.completed_chunks));
    frame_interior=manifest.positions(1,:)<315&manifest.positions(2,:)<315;
    interior=repmat(frame_interior(:),nf,1);
    candidate_mask=interior&w.code==1&w.unchanged;
    assert(sum(candidate_mask)==24433774&&sum(interior)==30005640);
    old=load(fullfile(inputDir,'masks.mat'),'common');
    assert(isequal(candidate_mask,old.common),'Common cohort differs from the original summary.');
    clear old w interior
    W=zeros(2,4,np);
    for p=1:np
        cable=WDM.param.M-manifest.positions(:,p);lengths=sqrt(sum(cable.^2,1));
        active=lengths>WDM.param.TOLL;W(:,active,p)=cable(:,active)./lengths(active);
    end
    savedGeometry=load(fullfile(inputDir,'geometry.mat'),'geometry');
    geometry_difference=max(abs(W(:)-savedGeometry.geometry.W(:)));
    assert(geometry_difference<=8*eps,'Reconstructed W differs from captured geometry.');
    clear savedGeometry
    files=cell(1,4);
    for k=1:4
        files{k}=matfile(fullfile(inputDir,[file_names{k} '.mat']));f=files{k};
        assert(isequal(size(f,'T'),[n 4])&&all(f.completed_chunks),'T is absent or incomplete.');
    end
    wecFile=matfile(fullfile(inputDir,'wec.mat'));
    pair_valid_mask=false(n,3);exclusion_reason_bits=zeros(n,3,'uint8');
    differences=nan(n,3);exclusive=zeros(3,8);overlap=zeros(3,7);
    readback=zeros(4,6);maxdiff=zeros(4,3);extremeRows=table();
    reason_legend={'DM output code or finite-tension failure','DM bound failure', ...
        'DM raw-force failure','method output code or finite-tension failure', ...
        'method bound failure','method raw-force failure','method native status failure'};
    chunk=250000;
    for first=1:chunk:n
        ix=(first:min(n,first+chunk-1))';local=find(candidate_mask(ix));
        if isempty(local),continue;end
        ids=ix(local);pid=mod(ids-1,np)+1;fid=floor((ids-1)/np)+1;
        F=manifest.requested_forces(:,fid)';
        fw=wecFile.Fw(ix,1:2);fw=fw(local,:);
        assert(all(sqrt(sum((fw-F).^2,2))<=1e-12),'Candidate force was not unchanged.');
        vx=reshape(W(1,:,pid),4,[])';vy=reshape(W(2,:,pid),4,[])';
        actual=cell(1,4);metrics=cell(1,4);
        for k=1:4
            f=files{k};T=f.T(ix,1:4);T=T(local,:);
            code=f.output_code(ix,1);code=code(local);
            native=ones(numel(ids),1);
            if native_applicable(k),native=f.diagnostics(ix,1);native=native(local);end
            q=physical_metrics(T,code,native,native_applicable(k),vx,vy,F,WDM);
            storedNorm=f.tension_norm_N(ix,1);storedNorm=storedNorm(local);
            storedBound=f.bound_violation_N(ix,1);storedBound=storedBound(local);
            storedRaw=f.force_error_raw_N(ix,1);storedRaw=storedRaw(local);
            acceptance=f.acceptance(ix,1);acceptance=acceptance(local);
            finiteT=all(isfinite(T),2);
            forward=128*eps*max(1,sqrt(sum([sum(abs(T.*vx),2),sum(abs(T.*vy),2)].^2,2)));
            differencesNorm=abs(q.norm-storedNorm);differencesBound=abs(q.bound-storedBound);
            differencesRaw=abs(q.error-storedRaw);
            normBad=finiteT&differencesNorm>8*eps(max(1,q.norm));
            boundBad=finiteT&differencesBound>8*eps(max(1,q.bound));
            rawBad=finiteT&differencesRaw>forward+8*eps(max(1,q.error));
            mismatch=q.valid~=logical(bitget(acceptance,6));
            readback(k,:)=readback(k,:)+[numel(ids),sum(normBad),sum(boundBad),sum(rawBad), ...
                sum(mismatch),sum(q.finite~=finiteT)];
            maxdiff(k,:)=max(maxdiff(k,:),[finite_max(differencesNorm),finite_max(differencesBound),finite_max(differencesRaw)]);
            assert(~any(normBad|boundBad|rawBad|mismatch),'Saved/independent physical readback differs.');
            actual{k}=T;metrics{k}=q;
        end
        d=metrics{1};
        for j=1:3
            q=metrics{j+1};reasons=zeros(numel(ids),1,'uint8');
            bad={~d.finite,d.finite&~d.bound_ok,d.finite&~d.force_ok, ...
                ~q.finite,q.finite&~q.bound_ok,q.finite&~q.force_ok,q.finite&~q.native_ok};
            remaining=true(numel(ids),1);
            for bit=1:7
                reasons=bitset(reasons,bit,bad{bit});overlap(j,bit)=overlap(j,bit)+sum(bad{bit});
                exclusive(j,bit)=exclusive(j,bit)+sum(remaining&bad{bit});remaining=remaining&~bad{bit};
            end
            valid=reasons==0;assert(isequal(valid,d.valid&q.valid));
            exclusive(j,8)=exclusive(j,8)+sum(valid);
            pair_valid_mask(ids,j)=valid;exclusion_reason_bits(ids,j)=reasons;
            delta=q.norm(valid)-d.norm(valid);differences(ids(valid),j)=delta;
            assert(all(isfinite(delta)),'A retained norm difference is nonfinite.');
            extremeRows=[extremeRows;extreme_rows(method_names{j},ids(valid),delta, ...
                actual{1}(valid,:),actual{j+1}(valid,:),d.error(valid),q.error(valid), ...
                q.bound(valid),q.native(valid))]; %#ok<AGROW>
        end
        if mod(floor((first-1)/chunk),20)==0
            fprintf('Saved-tension rows checked through %d/%d.\n',ix(end),n);
        end
    end
    assert(all(readback(:,1)==sum(candidate_mask))&&all(readback(:,2:end)==0,'all'));
    statistics=table();exclusions=table();finalExtremes=table();
    for j=1:3
        delta=differences(pair_valid_mask(:,j),j);N=numel(delta);
        qs=type5_quantiles(delta,[.025 .5 .975]);reference=quantile(delta,[.025 .5 .975]);
        assert(all(abs(qs-reference)<=8*eps(max(1,abs(qs)))),'Type-5 implementation differs from MATLAB quantile.');
        row=table(string(method_names{j}),string('DM'),string('norm(T_method,2)-norm(T_DM,2)'), ...
            sum(candidate_mask),N,sum(candidate_mask)-N,qs(2),qs(1),qs(3),min(delta),max(delta), ...
            sum(delta<0),sum(delta==0),sum(delta>0),sum(delta < -1e-8),sum(delta>1e-8), ...
            1e-9,1e-8,1e-12,native_applicable(j+1), ...
            string('Hyndman-Fan type 5; p_i=(i-0.5)/N; linear interpolation; endpoint clamping'), ...
            'VariableNames',{'method','reference','sign_convention','candidate_N','valid_pairs_N', ...
            'excluded_N','median_difference_N','q2_5_difference_N','q97_5_difference_N', ...
            'minimum_difference_N','maximum_difference_N','negative_count','exact_zero_count', ...
            'positive_count','below_minus_1e_8_N_count','above_plus_1e_8_N_count','bound_tolerance_N', ...
            'raw_force_tolerance_N','unchanged_cohort_tolerance_N','native_positive_flag_required', ...
            'quantile_estimator'});
        statistics=[statistics;row]; %#ok<AGROW>
        assert(N==exclusive(j,8)&&sum(exclusive(j,:))==sum(candidate_mask));
        for bit=1:7
            exclusions=[exclusions;table(string(method_names{j}),bit,string(reason_legend{bit}), ...
                overlap(j,bit),exclusive(j,bit), ...
                'VariableNames',{'method','reason_bit','reason','overlapping_failure_count','exclusive_first_failure_count'})]; %#ok<AGROW>
        end
        for direction={'lowest','highest'}
            r=extremeRows(strcmp(extremeRows.method,method_names{j})&strcmp(extremeRows.direction,direction{1}),:);
            if strcmp(direction{1},'lowest'),r=sortrows(r,{'difference_N','input_row'},{'ascend','ascend'});
            else,r=sortrows(r,{'difference_N','input_row'},{'descend','ascend'});end
            r=r(1:min(10,height(r)),:);r.pose_id=mod(r.input_row-1,np)+1;r.force_id=floor((r.input_row-1)/np)+1;
            r.x_mm=manifest.positions(1,r.pose_id)';r.y_mm=manifest.positions(2,r.pose_id)';
            r.requested_Fx_N=manifest.requested_forces(1,r.force_id)';r.requested_Fy_N=manifest.requested_forces(2,r.force_id)';
            finalExtremes=[finalExtremes;r]; %#ok<AGROW>
        end
    end
    clear differences extremeRows
    independent_readback=table(string(file_names'),readback(:,1),readback(:,2),readback(:,3), ...
        readback(:,4),readback(:,5),readback(:,6),maxdiff(:,1),maxdiff(:,2),maxdiff(:,3), ...
        'VariableNames',{'file','candidate_rows_checked','norm_mismatches','bound_mismatches', ...
        'raw_force_mismatches','validity_mask_mismatches','code_finiteness_mismatches', ...
        'max_norm_difference_from_saved_N','max_bound_difference_from_saved_N','max_raw_force_difference_from_saved_N'});
    mask_metadata=struct('ninputs',n,'nposes',np,'nforces',nf,'methods',{method_names}, ...
        'candidate_definition','frame interior AND checked WEC code 1 AND actual change <= 1e-12 N', ...
        'pair_definition','candidate AND both code 1/four finite actual tensions AND bounds <= 1e-9 N AND raw force error <= 1e-8 N AND QP native exitflag > 0', ...
        'reason_bits',{reason_legend},'rows_outside_candidate_have_reason_zero',true, ...
        'input_mapping','pose_id=mod(row-1,nposes)+1;force_id=floor((row-1)/nposes)+1', ...
        'bound_tolerance_N',1e-9,'force_tolerance_N',1e-8,'unchanged_tolerance_N',1e-12, ...
        'objective_difference_never_used_for_exclusion',true);
    save(fullfile(outputDir,'norm_difference_masks.mat'),'candidate_mask','pair_valid_mask', ...
        'exclusion_reason_bits','mask_metadata','statistics','exclusions','independent_readback','-v7.3');
    writetable(statistics,fullfile(outputDir,'norm_difference_statistics.csv'));
    writetable(exclusions,fullfile(outputDir,'norm_difference_exclusions.csv'));
    writetable(finalExtremes,fullfile(outputDir,'norm_difference_extremes.csv'));
    writetable(independent_readback,fullfile(outputDir,'norm_difference_readback.csv'));
    fprintf('Verifying immutable input hashes again.\n');
    after=check_files(inputDir,inventory,names);
    writetable(after,fullfile(outputDir,'norm_input_hashes_after.csv'));
    assert(all(after.passed)&&isequal(before.observed_sha256,after.observed_sha256));
    assert(strcmp(receipt.inventory_sha256,hash_file(inventoryFile)));
    assert(strcmp(sourceHash,hash_file(sourceFile))&&strcmp(sourceHash, ...
        hash_file(fullfile(sourceDir,'compute_norm_differences_saved.m'))));
    write_notes(outputDir,statistics,geometry_difference);
    receipt.status='passed';receipt.finished_at=datestr(now,30);
    receipt.original_inputs_unchanged=true;receipt.input_files_checked=height(before);
    receipt.source_file=sourceFile;receipt.source_sha256=sourceHash;
    receipt.candidate_N=sum(candidate_mask);receipt.retained_pairs=table2struct(statistics);
    receipt.physical_readback_passed=true;receipt.quantile_estimator=char(statistics.quantile_estimator(1));
    outputFiles={'norm_difference_statistics.csv','norm_difference_exclusions.csv', ...
        'norm_difference_masks.mat','norm_difference_extremes.csv','norm_difference_readback.csv','norm_difference_notes.txt'};
    outputs=struct('file',{},'bytes',{},'sha256',{});
    for k=1:numel(outputFiles)
        p=fullfile(outputDir,outputFiles{k});d=dir(p);
        outputs(end+1)=struct('file',outputFiles{k},'bytes',d.bytes,'sha256',hash_file(p)); %#ok<AGROW>
    end
    receipt.outputs=outputs;write_json(fullfile(outputDir,'norm_difference_receipt.json'),receipt);
    result=receipt;fprintf('Saved-data norm analysis PASSED.\n');disp(statistics(:,[1 4:9]));
catch exception
    receipt.status='failed';receipt.error_identifier=exception.identifier;receipt.error_message=exception.message;
    write_json(fullfile(outputDir,'norm_difference_receipt.json'),receipt);rethrow(exception)
end
end

function q=physical_metrics(T,code,native,needsNative,vx,vy,F,WDM)
q.finite=code==1&all(isfinite(T),2);
q.norm=sqrt(sum(T.^2,2));achieved=[sum(T.*vx,2),sum(T.*vy,2)];
q.error=sqrt(sum((achieved-F).^2,2));
q.bound=max([zeros(size(T,1),1),WDM.param.lim_inf(:)'-T,T-WDM.param.lim_sup(:)'],[],2);
q.bound(~all(isfinite(T),2))=NaN;
q.bound_ok=q.bound<=1e-9;q.force_ok=q.error<=1e-8;
q.native=native;q.native_ok=true(size(native));
if needsNative,q.native_ok=isfinite(native)&native>0;end
q.valid=q.finite&q.bound_ok&q.force_ok&q.native_ok;
end

function values=type5_quantiles(x,p)
x=sort(x(:));n=numel(x);assert(n>0&&all(isfinite(x)));
r=n*p+.5;lower=max(1,min(n,floor(r)));upper=max(1,min(n,ceil(r)));
weight=r-floor(r);values=(1-weight).*x(lower)'+weight.*x(upper)';
end

function rows=extreme_rows(method,ids,delta,dm,T,dmerror,error,bound,native)
rows=table();if isempty(delta),return;end
[~,order]=sort(delta);take=min(10,numel(delta));
for direction={'lowest','highest'}
    if strcmp(direction{1},'lowest'),ix=order(1:take);else,ix=order(end-take+1:end);end
    r=table(repmat(string(method),numel(ix),1),repmat(string(direction{1}),numel(ix),1), ...
        ids(ix),delta(ix),dm(ix,1),dm(ix,2),dm(ix,3),dm(ix,4), ...
        T(ix,1),T(ix,2),T(ix,3),T(ix,4),dmerror(ix),error(ix),bound(ix),native(ix), ...
        'VariableNames',{'method','direction','input_row','difference_N','DM_T1_N','DM_T2_N', ...
        'DM_T3_N','DM_T4_N','method_T1_N','method_T2_N','method_T3_N','method_T4_N', ...
        'DM_raw_force_error_N','method_raw_force_error_N','method_bound_violation_N','native_flag_or_not_applicable_one'});
    rows=[rows;r]; %#ok<AGROW>
end
end

function t=check_files(folder,inventory,names)
t=table();
for k=1:numel(names)
    name=names{k};ix=find(strcmp(inventory.relative_path,name));assert(isscalar(ix));
    p=fullfile(folder,name);d=dir(p);assert(isscalar(d),'Missing input: %s',p);
    actual=hash_file(p);expected=inventory.sha256(ix);passed=d.bytes==inventory.bytes(ix)&&strcmpi(actual,expected);
    t=[t;table(string(name),d.bytes,string(expected),string(actual),passed, ...
        'VariableNames',{'relative_path','bytes','expected_sha256','observed_sha256','passed'})]; %#ok<AGROW>
end
end
function write_notes(folder,statistics,geometryDifference)
fid=fopen(fullfile(folder,'norm_difference_notes.txt'),'w');assert(fid~=-1);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,['SAVED FINAL-RUN ACTUAL-TENSION NORM DIFFERENCES\n\n', ...
    'Delta = sqrt(sum(T_method.^2)) - sqrt(sum(T_DM.^2)), in N, using all four actual\n', ...
    'cable tensions. No normalization, excess-tension subtraction, squared objective\n', ...
    'or relative objective difference is used. Negative/zero values are retained.\n\n', ...
    'The start cohort contains 24,433,774 unchanged-force requests at interior poses.\n', ...
    'For each method separately, both outputs must be code 1/four finite tensions,\n', ...
    'within 3.8--25 N bounds to 1e-9 N, and reproduce the identical original requested\n', ...
    'force with Euclidean residual <= 1e-8 N. Bounded QP also needs positive native\n', ...
    'exitflag. DM/ICFM/VTDA have no native convergence flag in this saved capture;\n', ...
    'a finite algebraic output alone is not treated as physical success.\n\n', ...
    'W is reconstructed from captured anchor positions, pose and normalization\n', ...
    'policy. Actual WEC change <= 1e-12 N is rechecked for every candidate. Physical\n', ...
    'metrics are recomputed from T, checked against saved metric fields and saved\n', ...
    'exact-raw acceptance. Tolerance-threshold mismatches fail the audit; there is\n', ...
    'no hidden threshold adjustment or objective-based exclusion.\n\n', ...
    'Quantiles use explicit Hyndman-Fan type 5: order x(1)<=...<=x(N), h=N*p+0.5,\n', ...
    'linearly interpolate adjacent order statistics, clamping beyond endpoints.\n', ...
    'This is the default MATLAB quantile midpoint convention, checked numerically.\n', ...
    'The median and 2.5th/97.5th percentiles describe the full retained pairs.\n\n', ...
    'norm_difference_masks.mat stores candidate_mask(Nrows,1), pair_valid_mask\n', ...
    '(Nrows,3) in ICFM/QP/VTDA order and uint8 exclusion_reason_bits(Nrows,3).\n', ...
    'Bits may overlap. The exclusions CSV also reports mutually exclusive first\n', ...
    'failures in increasing bit order; these sum to candidate_N-retained_N.\n', ...
    'Rows outside the candidate mask have zero reason bits but are never valid.\n', ...
    'Saved row mapping: pose=mod(row-1,4096)+1; force=floor((row-1)/4096)+1.\n\n', ...
    'The extrema CSV retains 10 smallest and 10 largest retained differences per\n', ...
    'method with actual tensions, positions, forces and physical checks. Values\n', ...
    'within floating-point scale must not be interpreted as meaningful optimality\n', ...
    'violations. Larger discrepancies are retained for investigation, not removed.\n', ...
    'This analysis performs zero solver calls and zero new timing observations.\n', ...
    'All-attempt timing cohorts are unchanged and are not filtered to these pairs.\n\n']);
fprintf(fid,'Maximum reconstructed/captured W element difference: %.17g\n',geometryDifference);
for k=1:height(statistics)
    r=statistics(k,:);fprintf(fid,'%s: candidate %d; retained %d; excluded %d; median %.17gN; q2.5 %.17gN; q97.5 %.17gN.\n', ...
        r.method,r.candidate_N,r.valid_pairs_N,r.excluded_N,r.median_difference_N,r.q2_5_difference_N,r.q97_5_difference_N);
end
end
function value=finite_max(x)
x=x(isfinite(x));value=0;if ~isempty(x),value=max(x);end
end
function write_json(file,value)
fid=fopen(file,'w');assert(fid~=-1);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
function hash=hash_file(file)
fid=fopen(char(file),'rb');assert(fid~=-1);c=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while true,b=fread(fid,16*1024*1024,'*uint8');if isempty(b),break;end;md.update(b);end
hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
function p=canonical(p)
p=char(p);
if ~java.io.File(p).isAbsolute(),p=fullfile(pwd,p);end
p=char(java.io.File(p).getCanonicalPath());
end
function yes=same_or_child(parent,child)
parent=lower(canonical(parent));child=lower(canonical(child));yes=strcmp(parent,child)||startsWith(child,[parent filesep]);
end
function restore_environment(folder,oldpath)
cd(folder);path(oldpath);
end
