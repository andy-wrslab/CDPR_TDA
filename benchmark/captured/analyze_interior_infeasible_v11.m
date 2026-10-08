function report=analyze_interior_infeasible_v11(folder)
%ANALYZE_INTERIOR_INFEASIBLE_V11 Paired accuracy within the spatial interior.
% Read completed v11 observations only; no solver or timer is called. Keep
% both-successful and each one-sided failure denominator explicit. Quantiles
% use identical original rows for X-ACTA and WEC+DM. Write a separate CSV/MAT
% without altering capture, summary, or previously generated result tables.
% The 1e-8 N tie rule and table schema match T_complete_paths_matched.csv.

manifest_file=fullfile(folder,'manifest.mat');
s=load(manifest_file,'manifest');m=s.manifest;
assert(any(strcmp(m.status,{'raw_complete','complete'})),'Raw capture is incomplete.');
np=size(m.positions,2);n=np*size(m.requested_forces,2);
assert(n==m.ninputs && np==m.nposes,'Input grid does not match its manifest.');
audit_file=fullfile(folder,'independent_readback_audit.mat');
a=load(audit_file,'report');
assert(height(a.report)==7 && all(a.report.passed) && all(a.report.rows==n), ...
    'Complete physical readback is required before paired analysis.');
source_file=[mfilename('fullpath') '.m'];source_sha256=hash_file(source_file);
manifest_sha256=hash_file(manifest_file);started_at=datestr(now,30);
x=matfile(fullfile(folder,'xacta_raw.mat'));
d=matfile(fullfile(folder,'wecdm_raw.mat'));
w=matfile(fullfile(folder,'wec.mat'));
assert(all(x.completed_chunks)&&all(d.completed_chunks)&&all(w.completed_chunks), ...
    'A source series has pending chunks.');
assert(isequal(size(x,'force_error_raw_N'),[n 1]) && ...
    isequal(size(d,'force_error_raw_N'),[n 1]) && ...
    isequal(size(w,'raw_class'),[n 1]),'Source row counts differ.');

names={'frame_interior_raw_infeasible','frame_interior_wec_changed'};
definitions={'x < 315 and y < 315 and raw_class == 3', ...
    'x < 315 and y < 315 and wec.code == 1 and actual unchanged flag == false'};
chunk_size=100000;nchunks=ceil(n/chunk_size);errors=cell(2,nchunks);
counts=zeros(2,5); % inputs, both, only X-ACTA, only WEC+DM, neither
block=0;
for first=1:chunk_size:n
    block=block+1;ix=(first:min(n,first+chunk_size-1))';
    pose=mod(ix-1,np)+1;
    interior=(m.positions(1,pose)<315 & m.positions(2,pose)<315)';
    raw_class=w.raw_class(ix,1);allowed=w.code(ix,1)==1;
    changed=~w.unchanged(ix,1);
    xs=logical(bitget(x.acceptance(ix,1),4));
    ds=logical(bitget(d.acceptance(ix,1),4));
    xe=x.force_error_raw_N(ix,1);de=d.force_error_raw_N(ix,1);
    assert(all(isfinite(xe(xs)))&&all(isfinite(de(ds))), ...
        'Successful output has nonfinite raw-force error.');
    masks={interior&raw_class==3,interior&allowed&changed};
    for k=1:2
        mask=masks{k};pair=mask&xs&ds;
        counts(k,:)=counts(k,:)+[sum(mask),sum(pair),sum(mask&xs&~ds), ...
            sum(mask&~xs&ds),sum(mask&~xs&~ds)];
        errors{k,block}=[xe(pair),de(pair)];
    end
end

tie_tolerance_N=1e-8;rows=cell(2,1);
for k=1:2
    paired=vertcat(errors{k,:});errors(k,:)={[]};
    assert(size(paired,1)==counts(k,2),'Paired error count differs.');
    difference=paired(:,2)-paired(:,1);
    nx=sum(difference>tie_tolerance_N);nd=sum(difference < -tie_tolerance_N);
    nt=sum(abs(difference)<=tie_tolerance_N);
    assert(nx+nd+nt==counts(k,2) && sum(counts(k,2:5))==counts(k,1), ...
        'Paired success/accuracy categories are not exhaustive.');
    fractions=[NaN NaN NaN];if counts(k,2)>0,fractions=[nx nd nt]/counts(k,2);end
    rows{k}=table(names(k),definitions(k),counts(k,1),counts(k,2), ...
        counts(k,3),counts(k,4),counts(k,5),tie_tolerance_N, ...
        finite_quantile(paired(:,1),.5),finite_quantile(paired(:,1),.99), ...
        finite_quantile(paired(:,2),.5),finite_quantile(paired(:,2),.99), ...
        finite_quantile(difference,.5),nx,nd,nt,fractions(1),fractions(2),fractions(3), ...
        'VariableNames',{'subset','subset_definition','inputs','both_successful_bounded', ...
        'only_xacta_successful_bounded','only_wecdm_successful_bounded','neither_successful_bounded', ...
        'pair_tie_tolerance_N','xacta_median_raw_error_N','xacta_p99_raw_error_N', ...
        'wecdm_median_raw_error_N','wecdm_p99_raw_error_N', ...
        'median_paired_wecdm_minus_xacta_error_N','xacta_closer_count','wecdm_closer_count', ...
        'tied_count','xacta_closer_fraction','wecdm_closer_fraction','tied_fraction'});
end
report=vertcat(rows{:});
assert(strcmp(manifest_sha256,hash_file(manifest_file)) && ...
    strcmp(source_sha256,hash_file(source_file)),'Manifest or analysis source changed.');
provenance=struct('status','complete','source_file',source_file,'source_sha256',source_sha256, ...
    'source_manifest_sha256',manifest_sha256,'input_count',n, ...
    'physical_audit_sha256',hash_file(audit_file), ...
    'source_series',{{'wec.mat','xacta_raw.mat','wecdm_raw.mat'}}, ...
    'started_at',started_at,'completed_at',datestr(now,30), ...
    'policy','Same original rows; acceptance bit 4 for bounded success; errors relative to raw F; 1e-8 N paired tie rule; no timing changes');
writetable(report,fullfile(folder,'T_complete_paths_matched_interiors.csv'));
save(fullfile(folder,'T_complete_paths_matched_interiors.mat'),'report','provenance');
disp(report);
end

function q=finite_quantile(values,p)
values=values(isfinite(values));if isempty(values),q=NaN;else,q=quantile(values,p);end
end

function value=hash_file(filename)
fid=fopen(filename,'rb');assert(fid>=0,'Cannot read source evidence: %s',filename);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
digest=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    bytes=fread(fid,1048576,'*uint8');if ~isempty(bytes),digest.update(bytes);end
end
value=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end
