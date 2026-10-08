function result = Compare_Review_v11(config)
%COMPARE_REVIEW_V11 Fresh matched post-WEC and raw-force reviewer experiment.
%   config.profile='full' uses all 30,965,760 original grid inputs; 'smoke'
%   uses 192. Every retained series has one observation per attempted input.
%   Eight fresh series: WEC alone; five allocations after the identical
%   checked WEC force; X-ACTA on the raw force; directly timed WEC+DM.
%   Nonfinite/physically invalid shared WEC skips the five dependent calls.
%   The direct pipeline includes WEC and finite dispatch checks, attempts DM
%   on every finite WEC force, and is physically assessed after timing.
%   No zero substitution, timing replacement, re-timing, or best-of selection.
%   config.analyze=false captures raw data only; summarize_Review_v11(outdir)
%   can be run separately. Raw arrays require completed_chunks=true; pending
%   preallocated rows are undefined until their chunk has been written.
if nargin<1,config=struct();end
if ischar(config)||isstring(config),config=struct('profile',char(config));end
config=defaults(config);here=fileparts(mfilename('fullpath'));
if isempty(config.outdir),config.outdir=fullfile(here,'reviewer_xacta','results',['review_v11_' config.profile]);end
assert(~exist(fullfile(config.outdir,'manifest.mat'),'file'), ...
    'Compare_Review_v11:ExistingRun','Choose a fresh output folder: %s',config.outdir);
if ~exist(config.outdir,'dir'),mkdir(config.outdir);end
old_path=path;path_cleanup=onCleanup(@()path(old_path)); %#ok<NASGU>
addpath(here);
previous_diary=get(0,'Diary');previous_diary_file=get(0,'DiaryFile');
diary_cleanup=onCleanup(@()restore_diary(previous_diary,previous_diary_file)); %#ok<NASGU>
diary(fullfile(config.outdir,'run_log.txt'));
methods={'DM post-WEC','ICFM post-WEC','Bounded QP post-WEC','VTDA-L2 post-WEC', ...
    'ACTA post-WEC','X-ACTA raw','WEC + DM raw'};
files={'dm_postwec','icfm_postwec','qp_postwec','vtda_postwec','acta_postwec','xacta_raw','wecdm_raw'};
functions={'wt4_2024_minmax_v4','wt_pott_v2','wt_qp_bounded_v11', ...
    'wt_gouttefarde_v2','wt_acta','wt_xacta','WEC_v5'};
for k=1:numel(functions)
    assert(strcmpi(which(functions{k}),fullfile(here,[functions{k} '.m'])), ...
        'Compare_Review_v11:ShadowedSource','Unexpected source resolution: %s',functions{k});
end
WDM.param.M=[-315 315 315 -315;-315 -315 315 315];
WDM.param.TOLL=1e-14;WDM.param.TOLL_WEC=1e-9;
WDM.param.lim_inf=3.8*ones(4,1);WDM.param.lim_sup=25*ones(4,1);
WDM.verbose=false;WDM.acta=config.acta_options;WDM.xacta=config.xacta_options;
x=linspace(0,315,64);magnitude=linspace(0,10,21);
theta=linspace(0,2*pi,361);theta=theta(1:end-1);
[px,py]=meshgrid(x(config.position_indices),x(config.position_indices));
positions=[px(:)';py(:)'];
[a,r]=meshgrid(theta(config.angle_indices),magnitude(config.force_indices));a=a';r=r';
requested_forces=[r(:)'.*cos(a(:)');r(:)'.*sin(a(:)')];
nposes=size(positions,2);nforces=size(requested_forces,2);n=nposes*nforces;
chunk_size=min(config.chunk_size,max(1,ceil(n/numel(methods))));nchunks=ceil(n/chunk_size);
stream=RandStream('mt19937ar','Seed',config.seed);chunk_order=randperm(stream,nchunks);
geometry=make_geometry(WDM,positions);
save(fullfile(config.outdir,'geometry.mat'),'geometry','-v7.3');
manifest=struct('config',config,'methods',{methods},'files',{files}, ...
    'solver_functions',{functions},'positions',positions,'requested_forces',requested_forces, ...
    'ninputs',n,'nposes',nposes,'nforces',nforces,'chunk_size',chunk_size, ...
    'chunk_order',chunk_order,'method_order',zeros(nchunks,7), ...
    'started_at',datestr(now,30),'status','incomplete','matlab_version',version, ...
    'toolbox_versions',ver,'source_snapshot',snapshot_sources(here,config.outdir,functions), ...
    'input_order','pose varies fastest; then force magnitude and angle as original v8', ...
    'frame_interior','x<315 and y<315; zero axes are interior symmetry axes', ...
    'timing_policy','one raw wall-clock observation per attempted public call; observed maxima, not WCET', ...
    'shared_wec_policy','fresh WEC_v5; nonfinite or independently outside force skips all five allocations; no substitute', ...
    'pipeline_policy','fresh WEC_v5 plus finite check plus DM; physical checks after timer', ...
    'output_code_legend','1 finite; -1 exception; -2 malformed; -3 nonfinite; -4 shared-WEC skip; -5 invalid pipeline WEC; 0 pending', ...
    'acceptance_bits',{{'finite','bounds','convergence','successful','exact_target','exact_raw'}}, ...
    'wec_code_legend','1 finite and independently allowed; -1 exception; -2 malformed; -3 nonfinite; -4 finite outside', ...
    'domain_class_legend',{{'strictly_interior_beyond_roundoff','boundary_uncertain','outside_beyond_roundoff','invalid_or_degenerate'}}, ...
    'wec_change_flag_bits',{{'changed','magnitude_increased','magnitude_reduced','reversed','direction_changed','zero_to_nonzero','native_flag_disagrees'}}, ...
    'comparison_notes','WEC may increase magnitude or modify zero forces; native flag is not unchanged-force evidence', ...
    'qp_policy','Bounded QP retains both bounds; native positive exit flag required; no relaxed-upper-bound fallback', ...
    'tolerance_policy','Internal options recorded in config; ACTA absolute and X-ACTA scaled residuals differ; common physical checks use config force/bound tolerances', ...
    'pending_rows','Undefined unless corresponding completed_chunks is true; never summarize incomplete files');
save(fullfile(config.outdir,'manifest.mat'),'manifest','WDM','-v7.3');
pool=gcp('nocreate');
if isempty(pool),pool=parpool('Processes',config.workers);end
assert(~isa(pool,'parallel.ThreadPool')&&pool.NumWorkers==config.workers, ...
    'Compare_Review_v11:PoolMismatch','Use the requested number of process workers.');
worker_cleanup=onCleanup(@()restore_worker_warnings(pool)); %#ok<NASGU>
fetchOutputs(parfevalOnAll(pool,@worker_warning_state,1,'capture',here));
fetchOutputs(parfevalOnAll(pool,@verify_worker_solvers,1,here,functions));
pinning=pin_workers(pool,config.pin_pcores);
pin_cleanup=onCleanup(@()restore_affinity(pinning)); %#ok<NASGU>
manifest.hardware=struct('computer_name',getenv('COMPUTERNAME'), ...
    'processor_identifier',getenv('PROCESSOR_IDENTIFIER'),'pool_workers',pool.NumWorkers, ...
    'logical_processors',str2double(getenv('NUMBER_OF_PROCESSORS')), ...
    'computer',computer,'pinning',pinning);
manifest.warmup_calls_per_worker=fetchOutputs(parfevalOnAll(pool,@warm_worker,1,WDM));
initialize_files(config.outdir,files,n,nchunks);
manifest.measured_attempts=zeros(1,8);manifest.pipeline_shared_wec_mismatches=0;
save(fullfile(config.outdir,'manifest.mat'),'manifest','-append');
fprintf('Review v11: %s, %d poses x %d forces = %d inputs; %d chunks; one measured call/series.\n', ...
    config.profile,nposes,nforces,n,nchunks);
fprintf('Strict frame interior includes zero axes; boundary is x=315 or y=315.\n');
fprintf('Shared WEC invalid -> all five dependent calls skipped; X-ACTA still receives raw F.\n');
timer_run=tic;
for block=1:nchunks
    chunk=chunk_order(block);first=(chunk-1)*chunk_size+1;last=min(chunk*chunk_size,n);
    count=last-first+1;rows=(first:last)';
    pose_id=mod(rows-1,nposes)+1;force_id=floor((rows-1)/nposes)+1;
    P=positions(:,pose_id);F=requested_forces(:,force_id);
    permutation=randperm(stream,count);inverse=zeros(1,count);inverse(permutation)=1:count;
    Ptimed=P(:,permutation);Ftimed=F(:,permutation);
    fw_tmp=nan(count,2);wec_elapsed=nan(count,1);wec_code=zeros(count,1,'int8');native_scaled=-ones(count,1,'int8');
    parfor k=1:count
        [fw,scaled,dt,code]=timed_wec(WDM,Ptimed(:,k),Ftimed(:,k));
        fw_tmp(k,:)=fw';native_scaled(k)=scaled;wec_elapsed(k)=dt;wec_code(k)=code;
    end
    Fw=fw_tmp(inverse,:)';native_scaled=native_scaled(inverse);
    wec_elapsed=wec_elapsed(inverse);wec_code=wec_code(inverse);
    [raw_margin,raw_class,raw_nearest]=classify_forces(geometry,pose_id,F);
    [wec_margin,wec_class]=classify_forces(geometry,pose_id,Fw);
    wec_code(wec_code==1 & wec_class>=3)=int8(-4);
    allowed=wec_code==1;
    unchanged=all(isfinite(Fw),1)' & sqrt(sum((Fw-F).^2,1))'<=config.unchanged_tolerance;
    unchanged_exact=all(Fw==F,1)';
    changes=wec_changes(F,Fw,native_scaled,config.unchanged_tolerance);
    mf=matfile(fullfile(config.outdir,'wec.mat'),'Writable',true);
    mf.Fw(first:last,:)=Fw';mf.elapsed_s(first:last,1)=wec_elapsed;
    mf.native_scaled(first:last,1)=native_scaled;mf.code(first:last,1)=wec_code;
    mf.raw_margin_N(first:last,1)=raw_margin;mf.wec_margin_N(first:last,1)=wec_margin;
    mf.raw_class(first:last,1)=raw_class;mf.wec_class(first:last,1)=wec_class;
    mf.raw_nearest_force_N(first:last,1)=raw_nearest;
    mf.unchanged(first:last,1)=unchanged;mf.unchanged_exact(first:last,1)=unchanged_exact;
    mf.wec_change_flags(first:last,1)=changes;mf.completed_chunks(chunk,1)=true;
    manifest.measured_attempts(1)=manifest.measured_attempts(1)+count;
    order=circshift(1:7,[0 -mod(block-1,7)]);manifest.method_order(block,:)=order;
    for m=order
        T=nan(count,4);elapsed_s=nan(count,1);output_code=zeros(count,1,'int8');
        diagnostics=nan(count,20);solver_flags=-32768*ones(count,3,'int16');
        pipeline_fw=nan(count,2);
        if m<=5,target=Fw(:,permutation);attempt=allowed(permutation);
        else,target=Ftimed;attempt=true(count,1);end
        parfor k=1:count
            if ~attempt(k)
                output_code(k)=-4;
            else
                [out,dt,code,fw]=timed_method(m,WDM,Ptimed(:,k),target(:,k));
                elapsed_s(k)=dt;output_code(k)=code;
                if isstruct(out)&&isfield(out,'T')&&isnumeric(out.T)&&numel(out.T)==4&&isreal(out.T)
                    T(k,:)=reshape(out.T,1,4);
                end
                if m==3||m==5||m==6,diagnostics(k,:)=read_diagnostics(out);end
                solver_flags(k,:)=read_flags(out);
                if m==7,pipeline_fw(k,:)=fw';end
            end
        end
        T=T(inverse,:);elapsed_s=elapsed_s(inverse);output_code=output_code(inverse);
        diagnostics=diagnostics(inverse,:);solver_flags=solver_flags(inverse,:);
        if m==7
            pipeline_fw=pipeline_fw(inverse,:);
            equal=all((pipeline_fw==Fw') | (isnan(pipeline_fw)&isnan(Fw')),2);
            % A downstream DM exception still has the newly computed WEC F.
            % If WEC itself throws, both records are missing and agree.
            manifest.pipeline_shared_wec_mismatches=manifest.pipeline_shared_wec_mismatches+sum(~equal);
        end
        if m==6,physical_target=F;else,physical_target=Fw;end
        [raw_error,target_error,bound_violation,tension_norm,acceptance]= ...
            assess_outputs(T,output_code,diagnostics,m,geometry,pose_id,F,physical_target,WDM,config);
        mf=matfile(fullfile(config.outdir,[files{m} '.mat']),'Writable',true);
        mf.T(first:last,:)=T;mf.elapsed_s(first:last,1)=elapsed_s;
        mf.output_code(first:last,1)=output_code;mf.acceptance(first:last,1)=acceptance;
        mf.force_error_raw_N(first:last,1)=raw_error;mf.force_error_target_N(first:last,1)=target_error;
        mf.bound_violation_N(first:last,1)=bound_violation;mf.tension_norm_N(first:last,1)=tension_norm;
        mf.solver_flags(first:last,:)=solver_flags;
        if m==3||m==5||m==6,mf.diagnostics(first:last,:)=diagnostics;end
        if m==7,mf.pipeline_Fw(first:last,:)=pipeline_fw;end
        mf.completed_chunks(chunk,1)=true;
        manifest.measured_attempts(m+1)=manifest.measured_attempts(m+1)+sum(isfinite(elapsed_s));
    end
    manifest.completed_chunks=block;manifest.capture_seconds=toc(timer_run);
    save(fullfile(config.outdir,'manifest.mat'),'manifest','-append');
    fprintf('Chunk %d/%d complete, inputs %d:%d; WEC allowed %d/%d; %.1f min.\n', ...
        block,nchunks,first,last,sum(allowed),count,manifest.capture_seconds/60);
end
assert(manifest.pipeline_shared_wec_mismatches==0,'Compare_Review_v11:WECDifference', ...
    'Directly recomputed pipeline WEC differed from shared WEC; inspect raw pipeline_Fw.');
assert(verify_sources(manifest.source_snapshot),'Compare_Review_v11:SourceChanged', ...
    'Production source changed during measurement.');
manifest.status='raw_complete';manifest.finished_at=datestr(now,30);
manifest.capture_seconds=toc(timer_run);manifest.source_snapshot_unchanged=true;
save(fullfile(config.outdir,'manifest.mat'),'manifest','-append');
result=manifest;
fprintf('All raw capture complete: %.1f min; %s\n',manifest.capture_seconds/60,config.outdir);
if config.analyze,result=summarize_Review_v11(config.outdir);end
end

function c=defaults(c)
if ~isfield(c,'profile'),c.profile='smoke';end
switch lower(c.profile)
    case 'smoke',xi=[1 33 63 64];fi=[1 11 21];ai=[1 91 181 271];
    case 'representative',xi=unique(round(linspace(1,64,17)));fi=1:21;ai=1:5:360;
    case 'full',xi=1:64;fi=1:21;ai=1:360;
    otherwise,error('Unknown profile.');
end
d=struct('outdir','','workers',8,'pin_pcores',true,'chunk_size',100000, ...
    'seed',20261006,'position_indices',xi,'force_indices',fi,'angle_indices',ai, ...
    'force_tolerance',1e-8,'bound_tolerance',1e-9,'unchanged_tolerance',1e-12, ...
    'near_boundary_margin',1e-8,'acta_options',struct('tol',1e-10), ...
    'xacta_options',struct('tol',1e-9,'force_tol',1e-8),'analyze',false);
unknown=setdiff(fieldnames(c),[fieldnames(d);{'profile'}]);
assert(isempty(unknown),'Unknown option: %s',strjoin(unknown,', '));
names=fieldnames(d);for k=1:numel(names),if ~isfield(c,names{k}),c.(names{k})=d.(names{k});end,end
validateattributes(c.workers,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.chunk_size,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.seed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
validateattributes(c.position_indices,{'numeric'},{'vector','integer','>=',1,'<=',64});
validateattributes(c.force_indices,{'numeric'},{'vector','integer','>=',1,'<=',21});
validateattributes(c.angle_indices,{'numeric'},{'vector','integer','>=',1,'<=',360});
for f={'force_tolerance','bound_tolerance','unchanged_tolerance','near_boundary_margin'}
    validateattributes(c.(f{1}),{'numeric'},{'scalar','finite','positive'});
end
assert(numel(unique(c.position_indices))==numel(c.position_indices),'Duplicate positions.');
assert(numel(unique(c.force_indices))==numel(c.force_indices),'Duplicate magnitudes.');
assert(numel(unique(c.angle_indices))==numel(c.angle_indices),'Duplicate angles.');
end

function initialize_files(folder,files,n,nchunks)
completed_chunks=false(nchunks,1);diag_fields=diagnostic_fields();
flag_fields={'native_scale_dir','native_case','native_fallback'};
for m=1:7
    file=fullfile(folder,[files{m} '.mat']);save(file,'completed_chunks','flag_fields','-v7.3');
    f=matfile(file,'Writable',true);f.T(n,4)=NaN;f.elapsed_s(n,1)=NaN;
    f.output_code(n,1)=int8(0);f.acceptance(n,1)=uint8(0);
    f.force_error_raw_N(n,1)=NaN;f.force_error_target_N(n,1)=NaN;
    f.bound_violation_N(n,1)=NaN;f.tension_norm_N(n,1)=NaN;
    f.solver_flags(n,3)=int16(-32768);
    if m==3||m==5||m==6,f.diagnostics(n,20)=NaN;save(file,'diag_fields','-append');end
    if m==7,f.pipeline_Fw(n,2)=NaN;end
end
file=fullfile(folder,'wec.mat');save(file,'completed_chunks','-v7.3');f=matfile(file,'Writable',true);
f.Fw(n,2)=NaN;f.elapsed_s(n,1)=NaN;f.native_scaled(n,1)=int8(-1);f.code(n,1)=int8(0);
f.raw_margin_N(n,1)=NaN;f.wec_margin_N(n,1)=NaN;f.raw_nearest_force_N(n,1)=NaN;
f.raw_class(n,1)=uint8(0);f.wec_class(n,1)=uint8(0);
f.unchanged(n,1)=false;f.unchanged_exact(n,1)=false;f.wec_change_flags(n,1)=uint16(0);
end

function [fw,scaled,elapsed,code]=timed_wec(WDM,p,f)
fw=[NaN;NaN];scaled=int8(-1);code=int8(1);timer=[];
try
    timer=tic;[candidate,native]=WEC_v5(WDM,p,f);elapsed=toc(timer);
    if ~isnumeric(candidate)||~isreal(candidate)||numel(candidate)~=2
        code=int8(-2);
    else
        fw=candidate(:);
        if any(~isfinite(fw)),code=int8(-3);end
    end
    if isscalar(native)&&(isnumeric(native)||islogical(native)),scaled=int8(native);end
catch
    elapsed=toc(timer);code=int8(-1);
end
end

function [out,elapsed,code,fw]=timed_method(m,WDM,p,f)
out=struct();fw=[NaN;NaN];timer=[];code=int8(1);
try
    switch m
        case 1,timer=tic;out=wt4_2024_minmax_v4(WDM,p,f);elapsed=toc(timer);
        case 2,timer=tic;out=wt_pott_v2(WDM,p,f);elapsed=toc(timer);
        case 3,timer=tic;out=wt_qp_bounded_v11(WDM,p,f);elapsed=toc(timer);
        case 4,timer=tic;out=wt_gouttefarde_v2(WDM,p,f);elapsed=toc(timer);
        case 5,timer=tic;out=wt_acta(WDM,p,f);elapsed=toc(timer);
        case 6,timer=tic;out=wt_xacta(WDM,p,f);elapsed=toc(timer);
        case 7
            % The complete online path, including WEC's polygon construction
            % and finite-output dispatch, is inside this single interval.
            timer=tic;
            fw=WEC_v5(WDM,p,f);
            if isnumeric(fw)&&isreal(fw)&&numel(fw)==2&&all(isfinite(fw(:)))
                fw=fw(:);out=wt4_2024_minmax_v4(WDM,p,fw);
            else
                code=int8(-5);
                if ~isnumeric(fw)||~isreal(fw)||numel(fw)~=2,fw=[NaN;NaN];else,fw=fw(:);end
            end
            elapsed=toc(timer);
    end
catch
    elapsed=toc(timer);code=int8(-1);return
end
if code==-5,return;end
if ~isstruct(out)||~isfield(out,'T')||~isnumeric(out.T)||numel(out.T)~=4||~isreal(out.T)
    code=int8(-2);
elseif any(~isfinite(out.T(:))),code=int8(-3);end
end

function fields=diagnostic_fields()
fields={'exitflag','mu','iter','stage','rac_iters','eac_iters','resid', ...
    'raw_resid','primal_resid','stationarity_resid','rac_resid','eac_resid', ...
    'stationarity_tolerance','roundoff_floor','barrier_evals','rac_barrier_evals', ...
    'eac_barrier_evals','backtracks','rac_backtracks','eac_backtracks'};
end
function d=read_diagnostics(out)
fields=diagnostic_fields();d=nan(1,20);if ~isstruct(out),return;end
for k=1:20
    if ~isfield(out,fields{k}),continue;end
    v=out.(fields{k});
    if (isnumeric(v)||islogical(v))&&isscalar(v),d(k)=double(v);
    elseif k==4&&(ischar(v)||isstring(v))
        if any(strcmpi(v,{'acta','rac','regularized','rac_only'})),d(k)=1;
        elseif any(strcmpi(v,{'eac','extended','rac_eac'})),d(k)=2;end
    end
end
end
function d=read_flags(out)
d=-32768*ones(1,3,'int16');fields={'scale_dir','case','fallback'};
if ~isstruct(out),return;end
for k=1:3
    if isfield(out,fields{k})&&(isnumeric(out.(fields{k}))||islogical(out.(fields{k})))&&isscalar(out.(fields{k}))&&isfinite(out.(fields{k}))
        d(k)=int16(out.(fields{k}));
    end
end
end

function geometry=make_geometry(WDM,P)
% The exact four-generator zonotope has at most eight facets. Floating-
% point hulls can retain additional nearly collinear box images; preserve
% all 16 candidates rather than assuming exact symbolic simplification.
n=size(P,2);W=zeros(2,4,n);normals=zeros(2,16,n);offsets=inf(16,n);counts=zeros(n,1,'uint8');
vertices=cell(n,1);rank_ok=false(n,1);
bits=double(bitget(repmat(uint8(0:15)',1,4),repmat(1:4,16,1)))';
box=WDM.param.lim_inf+(WDM.param.lim_sup-WDM.param.lim_inf).*bits;
for p=1:n
    v=WDM.param.M-P(:,p);
    for k=1:4,d=norm(v(:,k));if d<=WDM.param.TOLL,v(:,k)=[0;0];else,v(:,k)=v(:,k)/d;end,end
    W(:,:,p)=v;if rank(v)<2,continue;end
    q=unique((v*box)','rows');h=convhull(q(:,1),q(:,2));q=q(h(1:end-1),:);
    e=q([2:end 1],:)-q;l=sqrt(sum(e.^2,2));N=[e(:,2)./l,-e(:,1)./l]';
    count=size(q,1);assert(count<=16,'Hull exceeds the number of box images.');
    normals(:,1:count,p)=N;offsets(1:count,p)=sum(N'.*q,2);counts(p)=uint8(count);
    vertices{p}=q;rank_ok(p)=true;
end
geometry=struct('W',W,'normals',normals,'offsets',offsets,'counts',counts,'vertices',{vertices},'rank_ok',rank_ok);
end

function [margin,classes,nearest]=classify_forces(g,pose_ids,F)
count=size(F,2);margin=nan(count,1);classes=4*ones(count,1,'uint8');nearest=nan(count,1);
nposes=numel(g.counts);first_pose=pose_ids(1);
for p=1:nposes
    idx=mod(p-first_pose,nposes)+1:nposes:count;
    if isempty(idx)||~g.rank_ok(p),continue;end
    valid=all(isfinite(F(:,idx)),1);idx=idx(valid);if isempty(idx),continue;end
    nf=double(g.counts(p));N=g.normals(:,1:nf,p);d=g.offsets(1:nf,p);force=F(:,idx);
    values=min(d-N'*force,[],1)';margin(idx)=values;
    tolerance=256*eps*max(1,max(abs(d))+sum(abs(force),1)');
    c=2*ones(numel(idx),1,'uint8');c(values>tolerance)=1;c(values < -tolerance)=3;classes(idx)=c;
    if nargout>=3
        nearest(idx)=0;outside=find(c==3);
        if ~isempty(outside)
            q=g.vertices{p};next=q([2:end 1],:);points=force(:,outside)';best=inf(numel(outside),1);
            for e=1:size(q,1)
                edge=next(e,:)-q(e,:);t=min(max(((points-q(e,:))*edge')/(edge*edge'),0),1);
                candidate=q(e,:)+t.*edge;best=min(best,sqrt(sum((points-candidate).^2,2)));
            end
            nearest(idx(outside))=best;
        end
    end
end
end

function flags=wec_changes(F,Fw,native,tol)
n=size(F,2);flags=zeros(n,1,'uint16');valid=all(isfinite(Fw),1)';
r=sqrt(sum(F.^2,1))';w=sqrt(sum(Fw.^2,1))';change=sqrt(sum((Fw-F).^2,1))'>tol;
increased=w-r>tol;reduced=r-w>tol;nonzero=r>tol;
dot=sum(F.*Fw,1)';reversed=nonzero & dot < -tol.*max(1,r);
perpendicular=nan(n,1);perpendicular(nonzero)=abs(F(1,nonzero).*Fw(2,nonzero)-F(2,nonzero).*Fw(1,nonzero))'./r(nonzero);
direction=nonzero & perpendicular>tol;zero_changed=~nonzero&w>tol;
conditions={change,increased,reduced,reversed,direction,zero_changed,(native==1)~=change};
for k=1:numel(conditions),flags=bitset(flags,k,valid&conditions{k});end
end

function [raw_error,target_error,violation,tension_norm,acceptance]=assess_outputs(T,code,d,m,g,pose_ids,F,target,WDM,c)
n=size(T,1);achieved=nan(n,2);nposes=numel(g.counts);first_pose=pose_ids(1);
for p=1:nposes
    idx=mod(p-first_pose,nposes)+1:nposes:n;if isempty(idx),continue;end
    achieved(idx,:)=(g.W(:,:,p)*T(idx,:)')';
end
raw_error=sqrt(sum((achieved-F').^2,2));target_error=sqrt(sum((achieved-target').^2,2));
violation=max([zeros(n,1),WDM.param.lim_inf'-T,T-WDM.param.lim_sup'],[],2);
violation(any(~isfinite(T),2))=NaN;tension_norm=sqrt(sum(T.^2,2));
finite=code==1;bounded=finite&violation<=c.bound_tolerance;converged=finite;
if m==3||m==5||m==6,converged=finite&d(:,1)>0;end
success=bounded&converged;exact_target=success&target_error<=c.force_tolerance;
exact_raw=success&raw_error<=c.force_tolerance;conditions={finite,bounded,converged,success,exact_target,exact_raw};
acceptance=zeros(n,1,'uint8');for k=1:6,acceptance=bitset(acceptance,k,conditions{k});end
end

function count=warm_worker(WDM)
P=[100 0 310;100 0 310];F=[5 0 10;5 0 10];count=0;
for rep=1:8
    for k=1:3
        [fw,~,~,code]=timed_wec(WDM,P(:,k),F(:,k));count=count+1;
        if code==1,for m=1:5,timed_method(m,WDM,P(:,k),fw);count=count+1;end,end
        timed_method(6,WDM,P(:,k),F(:,k));timed_method(7,WDM,P(:,k),F(:,k));count=count+2;
    end
end
end

function state = pin_workers(pool,enabled)
state=struct('requested',logical(enabled),'applied',false, ...
    'reason','disabled or hardware/pool does not match v8', ...
    'worker_pids',[],'worker_affinity_masks',[],'client_affinity_mask',[], ...
    'original_worker_masks',[],'original_client_mask',[]);
if ~(enabled && ispc && feature('numcores')==16 && ...
        str2double(getenv('NUMBER_OF_PROCESSORS'))==24 && pool.NumWorkers<=8)
    fprintf('P-core pinning skipped (disabled or hardware/pool does not match v8).\n');
    return
end
state.worker_pids=fetchOutputs(parfevalOnAll(pool,@feature,1,'getpid'));
state.original_worker_masks=zeros(size(state.worker_pids),'int64');
for k=1:numel(state.worker_pids)
    process=System.Diagnostics.Process.GetProcessById(int32(state.worker_pids(k)));
    state.original_worker_masks(k)=process.ProcessorAffinity.ToInt64();
end
client=System.Diagnostics.Process.GetCurrentProcess();
state.original_client_mask=client.ProcessorAffinity.ToInt64();
state.worker_affinity_masks=zeros(size(state.worker_pids),'int64');
try
    for k=1:numel(state.worker_pids)
        process=System.Diagnostics.Process.GetProcessById(int32(state.worker_pids(k)));
        process.ProcessorAffinity=System.IntPtr(int64(bitshift(3,2*(k-1))));
        state.worker_affinity_masks(k)=process.ProcessorAffinity.ToInt64();
    end
    client.ProcessorAffinity=System.IntPtr(int64(bitshift(255,16)));
    state.client_affinity_mask=client.ProcessorAffinity.ToInt64();
    state.applied=true; state.reason='dedicated P-core masks assigned; restored on exit';
catch err
    restore_affinity(state);
    rethrow(err)
end
end

function restore_affinity(state)
if isempty(state.original_client_mask), return; end
for k=1:numel(state.worker_pids)
    try
        process=System.Diagnostics.Process.GetProcessById(int32(state.worker_pids(k)));
        process.ProcessorAffinity=System.IntPtr(state.original_worker_masks(k));
    catch err
        warning('Compare_Review_v11:RestoreAffinity','Could not restore worker affinity: %s',err.message);
    end
end
try
    client=System.Diagnostics.Process.GetCurrentProcess();
    client.ProcessorAffinity=System.IntPtr(state.original_client_mask);
catch err
    warning('Compare_Review_v11:RestoreAffinity','Could not restore client affinity: %s',err.message);
end
end

function result = worker_warning_state(action,here)
persistent saved
result=true;
if strcmp(action,'capture')
    saved=struct('warnings',warning,'path',path);
    addpath(here); warning('off','all');
elseif ~isempty(saved)
    warning(saved.warnings); path(saved.path); saved=[];
end
end

function result=verify_worker_solvers(here,functions)
for k=1:numel(functions)
    assert(strcmpi(which(functions{k}),fullfile(here,[functions{k} '.m'])), ...
        'Compare_Review_v11:WorkerSolverPath','Worker has a shadowed solver: %s',functions{k});
end
result=true;
end

function restore_worker_warnings(pool)
try
    fetchOutputs(parfevalOnAll(pool,@worker_warning_state,1,'restore'));
catch err
    warning('Compare_Review_v11:RestoreWarnings','Could not restore worker warnings: %s',err.message);
end
end

function restore_diary(state,filename)
diary('off');
diary(filename);
if ~strcmpi(state,'on'),diary('off');end
end

function records=snapshot_sources(here,outdir,functions)
names=[{'Compare_Review_v11'},functions];
records=repmat(struct('original_path','','snapshot_path','','sha256',''),numel(names),1);
folder=fullfile(outdir,'source_snapshot');if ~exist(folder,'dir'),mkdir(folder);end
for k=1:numel(names)
    source=fullfile(here,[names{k} '.m']);target=fullfile(folder,[names{k} '.m']);copyfile(source,target);
    records(k)=struct('original_path',source,'snapshot_path',target,'sha256',sha256_file(target));
end
end
function yes=verify_sources(records)
yes=true;for k=1:numel(records)
    yes=yes&&exist(records(k).original_path,'file')~=0&&strcmp(sha256_file(records(k).original_path),records(k).sha256);
end
end
function value=sha256_file(filename)
fid=fopen(filename,'rb');assert(fid>=0,'Cannot read source: %s',filename);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
bytes=fread(fid,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
