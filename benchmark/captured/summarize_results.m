function summary = summarize_results(folder,make_plots)
%SUMMARIZE_RESULTS Analyze immutable one-pass reviewer capture.
%   No solver is called and no raw value is replaced. At most a few scalar
%   columns are loaded together; full diagnostic matrices are never loaded.
%   Timings include attempted failures; skipped calls have NaN timing and
%   explicit denominators. Component sums are labeled derived, not directly
%   timed complete paths. QP-relative L2 differences require matched forces
%   and physical acceptance; they are not independent optimality proofs.
if nargin<2,make_plots=true;end
analysis_started_at=datestr(now,30);analysis_timer=tic;
S=load(fullfile(folder,'manifest.mat'));manifest=S.manifest;WDM=S.WDM;
assert(any(strcmp(manifest.status,{'raw_complete','complete'})),'Capture is incomplete.');
n=manifest.ninputs;methods=manifest.methods;files=manifest.files;
w=load(fullfile(folder,'wec.mat'),'completed_chunks','code','raw_class','wec_class', ...
    'raw_margin_N','wec_margin_N','raw_nearest_force_N','unchanged','unchanged_exact','native_scaled','wec_change_flags');
assert(all(w.completed_chunks),'WEC has incomplete chunks.');
pose_edge=manifest.positions(1,:)==315|manifest.positions(2,:)==315;
frame_edge=repmat(pose_edge(:),manifest.nforces,1);
allowed=w.code==1;interior=~frame_edge;
common=interior&allowed&w.unchanged;
mask_names={'all','frame_interior','frame_edge','postwec_checked_interior', ...
    'common_unchanged','wec_changed','wec_invalid','raw_force_strict','raw_force_infeasible'};
masks={true(n,1),interior,frame_edge,interior&allowed,common, ...
    allowed&~w.unchanged,~allowed,w.raw_class==1,w.raw_class==3};
domain=table(mask_names',cellfun(@sum,masks)', ...
    'VariableNames',{'subset','inputs'});
common_exact=common;T_timing=table();T_quality=table();T_L2_matched=table();
T_iterative_domains=table();T_postwec_components=table();T_force_approximation=table();
qp=matfile(fullfile(folder,'qp_postwec.mat'));
qp_norm=qp.tension_norm_N(:,1);qp_acceptance=qp.acceptance(:,1);
qp_exact=logical(bitget(qp_acceptance,5));
wec_time=read_time(folder,'wec',n);
T_timing=[T_timing; timing_rows('WEC alone',wec_time,mask_names,masks)];
for m=1:7
    file=matfile(fullfile(folder,[files{m} '.mat']));
    assert(all(file.completed_chunks),'Incomplete method file: %s',files{m});
    assert(isequal(size(file,'elapsed_s'),[n,1]),'Method timing length differs.');
    elapsed=file.elapsed_s(:,1);code=file.output_code(:,1);accept=file.acceptance(:,1);
    T_timing=[T_timing;timing_rows(methods{m},elapsed,mask_names,masks)]; %#ok<AGROW>
    successful=logical(bitget(accept,4));exact_target=logical(bitget(accept,5));
    exact_raw=logical(bitget(accept,6));
    if m<=6,common_exact=common_exact&exact_target;end
    raw_error=file.force_error_raw_N(:,1);target_error=file.force_error_target_N(:,1);
    violation=file.bound_violation_N(:,1);tension_norm=file.tension_norm_N(:,1);
    T_quality=[T_quality;quality_rows(methods{m},code,accept,raw_error,target_error, ...
        violation,mask_names,masks)]; %#ok<AGROW>
    approximation_names={'all','frame_interior','frame_edge','raw_force_infeasible'};
    approximation_masks={true(n,1),interior,frame_edge,w.raw_class==3};
    for k=1:numel(approximation_names)
        mask=approximation_masks{k}&successful&isfinite(w.raw_nearest_force_N)&isfinite(raw_error);
        excess=raw_error(mask)-w.raw_nearest_force_N(mask);
        row=table(methods(m),approximation_names(k),sum(approximation_masks{k}),sum(mask), ...
            finite_quantile(raw_error(mask),.5),finite_quantile(w.raw_nearest_force_N(mask),.5), ...
            finite_quantile(excess,.5),finite_quantile(excess,.99),finite_min(excess),finite_max(excess), ...
            sum(excess>manifest.config.force_tolerance),sum(excess < -manifest.config.force_tolerance), ...
            'VariableNames',{'method','subset','inputs','successful_with_reference', ...
            'median_achieved_raw_error_N','median_independent_nearest_distance_N', ...
            'median_excess_over_nearest_N','p99_excess_over_nearest_N','minimum_excess_over_nearest_N', ...
            'maximum_excess_over_nearest_N','excess_above_force_tolerance','below_nearest_beyond_force_tolerance'});
        T_force_approximation=[T_force_approximation;row]; %#ok<AGROW>
    end
    if m<=5
        mask=interior&allowed;
        component=elapsed(mask)+wec_time(mask);
        a=statistics_us(elapsed(mask));b=statistics_us(wec_time(mask));c=statistics_us(component);
        row=table(methods(m),sum(mask),sum(isfinite(component)), ...
            b(3),a(3),c(3),c(5),c(7), ...
            'VariableNames',{'method','inputs','finite_component_sums','shared_wec_median_us', ...
            'allocation_median_us','derived_sum_median_us','derived_sum_p99_us','derived_sum_observed_max_us'});
        T_postwec_components=[T_postwec_components;row]; %#ok<AGROW>
    end
    % The raw X-ACTA target matches the shared WEC target only on unchanged
    % requests. Every other allocator uses that same target directly.
    matched=interior&allowed&qp_exact&exact_target;
    if m==6,matched=matched&w.unchanged&exact_raw;end
    delta=tension_norm(matched)-qp_norm(matched);
    objective_gap=.5*(tension_norm(matched).^2-qp_norm(matched).^2);
    row=table(methods(m),sum(matched),finite_quantile(delta,.5),finite_quantile(delta,.99), ...
        finite_min(delta),finite_max(delta),finite_quantile(objective_gap,.5), ...
        finite_max(abs(objective_gap)),sum(delta < -1e-8), ...
        'VariableNames',{'method','matched_physically_accepted_inputs','median_norm_difference_N', ...
        'p99_norm_difference_N','min_norm_difference_N','max_norm_difference_N', ...
        'median_half_norm_squared_difference_N2','max_absolute_half_norm_squared_difference_N2', ...
        'outputs_below_QP_norm_by_more_than_1e_8_N'});
    T_L2_matched=[T_L2_matched;row]; %#ok<AGROW>
    if any(m==[3 5 6])
        native_flag=file.diagnostics(:,1);iterations=file.diagnostics(:,3);
        stage=file.diagnostics(:,4);residual=file.diagnostics(:,7);
        if m==6,force_class=w.raw_class;margin=w.raw_margin_N;
        else,force_class=w.wec_class;margin=w.wec_margin_N;end
        iterative_names={'all','frame_interior','frame_edge','strict_away_from_boundary', ...
            'strict_thin_interior','boundary_uncertain','force_outside','common_unchanged'};
        iterative_masks={true(n,1),interior,frame_edge, ...
            interior&force_class==1&margin>manifest.config.near_boundary_margin, ...
            interior&force_class==1&margin<=manifest.config.near_boundary_margin, ...
            interior&force_class==2,interior&force_class==3,common};
        for k=1:numel(iterative_names)
            mask=iterative_masks{k};attempt=mask&isfinite(elapsed);
            good=attempt&successful;
            row=table(methods(m),iterative_names(k),sum(mask),sum(attempt),sum(attempt&native_flag>0), ...
                sum(good),sum(mask&exact_target),sum(attempt&stage==1),sum(attempt&stage==2), ...
                finite_mean(iterations(attempt)),finite_quantile(iterations(attempt),.5),finite_quantile(iterations(attempt),.99), ...
                finite_max(iterations(attempt)),finite_quantile(residual(attempt),.5), ...
                finite_max(residual(good)), ...
                'VariableNames',{'method','subset','inputs','attempts','native_converged', ...
                'successful_bounded','exact_target','stage1','stage2','mean_iterations','median_iterations', ...
                'p99_iterations','max_iterations','median_native_residual','max_native_residual_successful'});
            T_iterative_domains=[T_iterative_domains;row]; %#ok<AGROW>
        end
        clear native_flag iterations stage residual
    end
    clear elapsed code accept successful exact_target exact_raw raw_error target_error violation tension_norm
    fprintf('Summary: %s complete.\n',methods{m});
end
% Publish exact denominators separately; the main common timing comparison
% remains on common unchanged requests even when an allocator fails.
domain=[domain;table({'common_exact_outputs'},sum(common_exact), ...
    'VariableNames',{'subset','inputs'})];
save(fullfile(folder,'masks.mat'),'common','common_exact','frame_edge','allowed','-v7.3');
T_wec_behavior=wec_behavior(w,interior);
T_pipeline_parity=pipeline_parity(folder,n,allowed);
T_complete_paths_matched=complete_paths_matched(folder,n,interior,w);
summary=struct('manifest',manifest,'T_domain',domain,'T_timing',T_timing, ...
    'T_quality',T_quality,'T_postwec_components',T_postwec_components, ...
    'T_L2_matched',T_L2_matched,'T_iterative_domains',T_iterative_domains, ...
    'T_wec_behavior',T_wec_behavior,'T_pipeline_parity',T_pipeline_parity, ...
    'T_force_approximation',T_force_approximation,'T_complete_paths_matched',T_complete_paths_matched, ...
    'analysis_source',analysis_sources(folder),'analyzed_at',datestr(now,30));
names={'T_domain','T_timing','T_quality','T_postwec_components','T_L2_matched', ...
    'T_iterative_domains','T_wec_behavior','T_pipeline_parity','T_force_approximation','T_complete_paths_matched'};
for k=1:numel(names),writetable(summary.(names{k}),fullfile(folder,[names{k} '.csv']));end
save(fullfile(folder,'summary.mat'),'summary','-v7.3');
write_report(folder,summary);
if make_plots,plot_results(folder,summary);end
% Keep the original measured manifest immutable. This separate receipt is
% written last, after every requested analysis artifact succeeds.
analysis_receipt=struct('status','analysis_complete','capture_status',manifest.status, ...
    'started_at',analysis_started_at,'completed_at',datestr(now,30), ...
    'elapsed_s',toc(analysis_timer),'plots_requested',logical(make_plots), ...
    'analysis_source',summary.analysis_source,'ninputs',n, ...
    'paired_complete_path_table','T_complete_paths_matched.csv');
save(fullfile(folder,'analysis_receipt.mat'),'analysis_receipt');
end

function t=read_time(folder,name,n)
f=matfile(fullfile(folder,[name '.mat']));assert(all(f.completed_chunks),'Incomplete series.');
assert(isequal(size(f,'elapsed_s'),[n,1]),'Time shape mismatch.');t=f.elapsed_s(:,1);
end

function rows=timing_rows(method,elapsed,names,masks)
rows=table();
for k=1:numel(names)
    t=elapsed(masks{k});s=statistics_us(t);
    row=table({method},names(k),numel(t),sum(isfinite(t)),sum(isnan(t)), ...
        s(1),s(2),s(3),s(4),s(5),s(6),s(7),s(8), ...
        'VariableNames',{'method','subset','inputs','timed_attempts','skipped_or_missing', ...
        'mean_us','q25_us','median_us','q75_us','p99_us','p99_9_us','observed_max_us','observed_min_us'});
    rows=[rows;row]; %#ok<AGROW>
end
end

function s=statistics_us(t)
t=t(isfinite(t))*1e6;
if isempty(t),s=nan(1,8);return;end
q=quantile(t,[.25 .5 .75 .99 .999]);s=[mean(t),q,max(t),min(t)];
end

function rows=quality_rows(method,code,accept,raw,target,violation,names,masks)
rows=table();success=logical(bitget(accept,4));
for k=1:numel(names)
    mask=masks{k};good=mask&success;finite=mask&code==1;
    row=table({method},names(k),sum(mask),sum(mask&code~=-4), ...
        sum(mask&code==-4),sum(mask&code==-5),sum(mask&code==-1), ...
        sum(mask&code==-2),sum(mask&code==-3),sum(finite), ...
        sum(mask&logical(bitget(accept,2))),sum(mask&logical(bitget(accept,3))), ...
        sum(good),sum(mask&logical(bitget(accept,5))),sum(mask&logical(bitget(accept,6))), ...
        finite_quantile(raw(good),.5),finite_quantile(raw(good),.99),finite_max(raw(good)), ...
        finite_quantile(target(good),.5),finite_max(target(good)),finite_max(violation(finite)), ...
        'VariableNames',{'method','subset','inputs','attempts','skipped_shared_wec', ...
        'pipeline_wec_failures','exceptions','malformed','nonfinite','finite_outputs', ...
        'bounded','native_converged_or_finite_algebraic','successful_bounded','exact_target','exact_raw', ...
        'successful_median_raw_error_N','successful_p99_raw_error_N','successful_max_raw_error_N', ...
        'successful_median_target_error_N','successful_max_target_error_N','max_bound_violation_all_finite_N'});
    rows=[rows;row]; %#ok<AGROW>
end
end

function rows=wec_behavior(w,interior)
names={'all','frame_interior','frame_edge'};masks={true(size(interior)),interior,~interior};rows=table();
for k=1:numel(names)
    mask=masks{k};row=table(names(k),sum(mask),sum(mask&w.code==1), ...
        sum(mask&w.code==-1),sum(mask&w.code==-2),sum(mask&w.code==-3),sum(mask&w.code==-4), ...
        sum(mask&w.native_scaled==1),sum(mask&w.unchanged),sum(mask&w.unchanged_exact), ...
        sum(mask&logical(bitget(w.wec_change_flags,2))),sum(mask&logical(bitget(w.wec_change_flags,3))), ...
        sum(mask&logical(bitget(w.wec_change_flags,4))),sum(mask&logical(bitget(w.wec_change_flags,5))), ...
        sum(mask&logical(bitget(w.wec_change_flags,6))),sum(mask&logical(bitget(w.wec_change_flags,7))), ...
        sum(mask&w.wec_class==1),sum(mask&w.wec_class==2),sum(mask&w.wec_class==3), ...
        finite_min(w.wec_margin_N(mask&w.code==1)), ...
        'VariableNames',{'subset','inputs','checked_allowed','exceptions','malformed','nonfinite', ...
        'finite_outside_polygon','native_scaled_true','actual_unchanged','exactly_unchanged', ...
        'magnitude_increased','magnitude_reduced','direction_reversed','direction_changed', ...
        'zero_to_nonzero','native_flag_disagrees_with_actual_change','strict_force_interior', ...
        'force_boundary_uncertain','force_outside','minimum_allowed_polygon_margin_N'});
    rows=[rows;row]; %#ok<AGROW>
end
end

function t=pipeline_parity(folder,n,allowed)
a=matfile(fullfile(folder,'dm_postwec.mat'));b=matfile(fullfile(folder,'wecdm_raw.mat'));
counts=zeros(1,4);maximum=0;
for first=1:100000:n
    ix=first:min(n,first+99999);mask=allowed(ix);A=a.T(ix,1:4);B=b.T(ix,1:4);
    ca=a.output_code(ix,1);cb=b.output_code(ix,1);
    finite=mask&all(isfinite(A),2)&all(isfinite(B),2);delta=max(abs(A(finite,:)-B(finite,:)),[],2);
    equal=all((A==B)|(isnan(A)&isnan(B)),2);
    counts=counts+[sum(mask),sum(mask&ca~=cb),sum(mask&~equal),sum(delta>1e-9)];
    if ~isempty(delta),maximum=max(maximum,max(delta));end
end
t=table(counts(1),counts(2),counts(3),counts(4),maximum, ...
    'VariableNames',{'checked_shared_inputs','output_code_changes','non_bit_identical_tensions', ...
    'finite_tension_changes_over_1e_9_N','maximum_finite_tension_change_N'});
end

function rows=complete_paths_matched(folder,n,interior,w)
% The pairwise accuracy denominator is identical for both methods. Keep
% non-overlapping success counts for the underlying subset so excluded
% failures remain visible. This does not filter any timing distribution.
x=matfile(fullfile(folder,'xacta_raw.mat'));d=matfile(fullfile(folder,'wecdm_raw.mat'));
xs=logical(bitget(x.acceptance(:,1),4));ds=logical(bitget(d.acceptance(:,1),4));
xe=x.force_error_raw_N(:,1);de=d.force_error_raw_N(:,1);
assert(numel(xs)==n&&numel(ds)==n&&numel(xe)==n&&numel(de)==n,'Pairwise input lengths differ.');
assert(all(isfinite(xe(xs)))&&all(isfinite(de(ds))),'Successful output has nonfinite raw-force error.');
names={'all','frame_interior','raw_infeasible','wec_changed'};
definitions={'all original rows','x < 315 and y < 315', ...
    'raw_class == 3; includes frame edge','wec.code == 1 and actual unchanged flag == false; includes frame edge'};
masks={true(n,1),interior,w.raw_class==3,w.code==1&~w.unchanged};
tie_tolerance_N=1e-8;rows=table();
for k=1:numel(names)
    mask=masks{k};pair=mask&xs&ds;count=sum(pair);
    paired_difference=de(pair)-xe(pair);
    nx=sum(paired_difference>tie_tolerance_N);nd=sum(paired_difference < -tie_tolerance_N);
    nt=sum(abs(paired_difference)<=tie_tolerance_N);
    assert(nx+nd+nt==count,'Paired error classification is not exhaustive.');
    fractions=[NaN NaN NaN];if count>0,fractions=[nx nd nt]/count;end
    row=table(names(k),definitions(k),sum(mask),count,sum(mask&xs&~ds), ...
        sum(mask&~xs&ds),sum(mask&~xs&~ds),tie_tolerance_N, ...
        finite_quantile(xe(pair),.5),finite_quantile(xe(pair),.99), ...
        finite_quantile(de(pair),.5),finite_quantile(de(pair),.99), ...
        finite_quantile(paired_difference,.5),nx,nd,nt,fractions(1),fractions(2),fractions(3), ...
        'VariableNames',{'subset','subset_definition','inputs','both_successful_bounded', ...
        'only_xacta_successful_bounded','only_wecdm_successful_bounded','neither_successful_bounded', ...
        'pair_tie_tolerance_N','xacta_median_raw_error_N','xacta_p99_raw_error_N', ...
        'wecdm_median_raw_error_N','wecdm_p99_raw_error_N', ...
        'median_paired_wecdm_minus_xacta_error_N','xacta_closer_count','wecdm_closer_count', ...
        'tied_count','xacta_closer_fraction','wecdm_closer_fraction','tied_fraction'});
    assert(row.both_successful_bounded+row.only_xacta_successful_bounded+ ...
        row.only_wecdm_successful_bounded+row.neither_successful_bounded==row.inputs, ...
        'Pairwise success categories do not cover the input subset.');
    rows=[rows;row]; %#ok<AGROW>
end
end

function q=finite_quantile(x,p)
x=x(isfinite(x));if isempty(x),q=NaN;else,q=quantile(x,p);end
end
function y=finite_mean(x)
x=x(isfinite(x));if isempty(x),y=NaN;else,y=mean(x);end
end
function y=finite_min(x)
x=x(isfinite(x));if isempty(x),y=NaN;else,y=min(x);end
end
function y=finite_max(x)
x=x(isfinite(x));if isempty(x),y=NaN;else,y=max(x);end
end

function records=analysis_sources(folder)
here=fileparts(mfilename('fullpath'));names={'summarize_results.m'};
destination=fullfile(folder,'analysis_source_snapshot');if ~exist(destination,'dir'),mkdir(destination);end
records=repmat(struct('path','','snapshot_path','','sha256',''),numel(names),1);
for k=1:numel(names)
    file=fullfile(here,names{k});records(k).path=file;
    assert(exist(file,'file')==2,'Analysis source is missing: %s',file);
    records(k).snapshot_path=fullfile(destination,names{k});copyfile(file,records(k).snapshot_path);
    fid=fopen(file,'rb');cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
    data=fread(fid,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(data);
    records(k).sha256=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));clear cleanup
end
end

function write_report(folder,s)
fid=fopen(fullfile(folder,'REVIEW_RESULTS.md'),'w');assert(fid>=0,'Cannot write report.');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Fresh single-pass reviewer comparison\n\n');
fprintf(fid,'%d inputs: %d poses and %d raw force requests per pose. All WEC and solver observations are fresh in the same capture session. Every attempted call retains its sole raw wall-clock observation, including failures and scheduling stalls. No fastest-of-five selection or replacement was made.\n\n',s.manifest.ninputs,s.manifest.nposes,s.manifest.nforces);
fprintf(fid,'The five allocation methods receive the same freshly computed, independently checked WEC force. Their timers include public-call geometry, initialization, solution and native diagnostics. Shared WEC is separately timed. Nonfinite or physically outside shared WEC output skips all five dependent calls, explicitly, without substituting zero.\n\n');
fprintf(fid,'X-ACTA receives raw requests. Its counterpart is the directly timed WEC+DM path, including WEC polygon construction and finite-output dispatch. Every finite WEC output is attempted by that complete path; independent physical checks follow the timer. The main six-method comparison uses frame-interior inputs whose checked WEC force is actually unchanged within %.3g N, retaining solver failures in its timing distribution.\n\n',s.manifest.config.unchanged_tolerance);
fprintf(fid,'| Subset | Inputs |\n|---|---:|\n');
for k=1:height(s.T_domain),fprintf(fid,'| %s | %d |\n',s.T_domain.subset{k},s.T_domain.inputs(k));end
fprintf(fid,'\n## Common unchanged-force timing\n\n| Method | Attempts | Median us | p99 us | Observed max us |\n|---|---:|---:|---:|---:|\n');
t=s.T_timing(strcmp(s.T_timing.subset,'common_unchanged') & ~strcmp(s.T_timing.method,'WEC alone') & ~strcmp(s.T_timing.method,'WEC + DM raw'),:);
for k=1:height(t),fprintf(fid,'| %s | %d | %.3f | %.3f | %.3f |\n',t.method{k},t.timed_attempts(k),t.median_us(k),t.p99_us(k),t.observed_max_us(k));end
fprintf(fid,'\n## Complete paths on all raw requests\n\n| Method | Attempts | Median us | p99 us | Observed max us |\n|---|---:|---:|---:|---:|\n');
t=s.T_timing(strcmp(s.T_timing.subset,'all') & ismember(s.T_timing.method,{'X-ACTA raw','WEC + DM raw'}),:);
for k=1:height(t),fprintf(fid,'| %s | %d | %.3f | %.3f | %.3f |\n',t.method{k},t.timed_attempts(k),t.median_us(k),t.p99_us(k),t.observed_max_us(k));end
fprintf(fid,'\nObserved maxima are measurements on a general-purpose operating system, not worst-case execution guarantees. T_postwec_components.csv shows sums of separately observed WEC and allocation times, explicitly derived component totals; those are not directly measured complete-path times. No timer includes offline reference polygons, independent output checks, input-grid construction, source hashing, storage, or plotting.\n\n');
fprintf(fid,'## Paired complete-path force accuracy\n\n');
fprintf(fid,'T_complete_paths_matched.csv compares raw-force error only on identical original rows where both complete paths return successful bounded outputs. It separately counts only-X-ACTA, only-WEC+DM and neither-successful outcomes within every underlying domain, so accuracy conditioning does not hide the failure denominators. A positive paired WEC+DM-minus-X-ACTA difference favors X-ACTA. Closer/tied classifications use a declared 1e-8 N threshold. This pairwise filter does not change the all-attempt timing comparison.\n\n');
fprintf(fid,'| Subset | Both successful | Only X-ACTA | Only WEC+DM | Neither | X-ACTA median error N | WEC+DM median error N | Median paired WEC+DM minus X-ACTA N |\n|---|---:|---:|---:|---:|---:|---:|---:|\n');
t=s.T_complete_paths_matched;
for k=1:height(t)
    fprintf(fid,'| %s | %d | %d | %d | %d | %.6g | %.6g | %.6g |\n', ...
        t.subset{k},t.both_successful_bounded(k),t.only_xacta_successful_bounded(k), ...
        t.only_wecdm_successful_bounded(k),t.neither_successful_bounded(k), ...
        t.xacta_median_raw_error_N(k),t.wecdm_median_raw_error_N(k),t.median_paired_wecdm_minus_xacta_error_N(k));
end
fprintf(fid,'\nThe raw_infeasible and wec_changed rows include frame-edge requests, while frame_interior excludes them. WEC-changed means independently allowed shared WEC output that differs from the raw request beyond the unchanged tolerance; it does not mean direction-preserving reduction only. Exact subset definitions and p99 errors are retained in the CSV. Per-method successful-error statistics elsewhere have their own denominators and cannot alone establish a pairwise accuracy advantage.\n\n');
fprintf(fid,'## Interpretation and acceptance\n\n');
fprintf(fid,'The spatial interior is x<315 mm and y<315 mm on the sampled nonnegative quadrant. Zero-coordinate symmetry axes are interior. Requests on x=315 or y=315 are reported separately. Wrench-domain classification uses independent 16-vertex tension-box image polygons and a scale-aware roundoff allowance. Actual signed margins are retained, so a tiny positive WEC margin is distinguishable from an uncertain boundary.\n\n');
fprintf(fid,'A finite tension is independently checked against bounds (%.3g N tolerance) and its requested target (%.3g N). Bounded QP, ACTA and X-ACTA also need their native positive exit flag. Bounded QP never drops the upper bounds. ACTA internal tol is %.3g; X-ACTA internal tol is %.3g with its returned roundoff allowance. These different residuals are not equated; common physical checks provide the comparison.\n\n', ...
    s.manifest.config.bound_tolerance,s.manifest.config.force_tolerance,s.manifest.config.acta_options.tol,s.manifest.config.xacta_options.tol);
fprintf(fid,'WEC is retained as implemented: it can increase force magnitude, modify zero force, or fail. Its native scaled flag remains true for some unchanged zero requests. T_wec_behavior.csv distinguishes these outcomes using the actual returned force. No statement that all WEC modifications are downward scaling is justified.\n\n');
fprintf(fid,'L2 differences in T_L2_matched.csv use only physically accepted outputs on the same force. X-ACTA is compared with post-WEC QP only where that force remains unchanged. ACTA/X-ACTA optimize centering objectives and need not minimize tension L2. QP-relative differences are comparative evidence, not independent optimality certificates; the separate active-set oracle validation provides that check on its documented sample.\n\n');
fprintf(fid,'T_force_approximation.csv compares each successful bounded output with the independent Euclidean distance from the raw request to the tension-box image polygon. Excess error is descriptive: WEC ray adjustment and analytic centering need not attain the nearest-point solution. Native iteration summaries distinguish the frame interior from its edge; the strict/thin/uncertain force-domain rows are restricted to the frame interior.\n\n');
fprintf(fid,'Iteration-domain bins refer to the force supplied to each method: shared post-WEC force for QP/ACTA and raw force for X-ACTA. These force-domain rows describe each method''s workload, not matched-request iteration comparisons. The common_unchanged row uses the shared matched input subset.\n\n');
fprintf(fid,'Timing and quality tables give all input, attempted, skipped, failure and successful-output denominators. Error quantiles labeled successful use bounded/converged outputs; failed last iterates remain in raw T and the failure counts. Data files preserve diagnostic residuals, iteration counts, native flags, and source versions.\n\n');
fprintf(fid,'Capture: %s to %s; MATLAB %s; %.1f min. Source hashes are in manifest.source_snapshot and analysis hashes in summary.analysis_source.\n',s.manifest.started_at,s.manifest.finished_at,s.manifest.matlab_version,s.manifest.capture_seconds/60);
end
