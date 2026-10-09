function WDM = wt_qp_bounded(WDM,P,F)
% WT_QP_BOUNDED Bounded minimum-norm cable-tension distribution.
%   OUT = wt_qp_bounded(WDM,P,F) solves min 0.5*T'*T subject to W*T=F
%   and lim_inf<=T<=lim_sup using MATLAB quadprog. Both bounds are retained
%   on every solve; there is no upper-bound-removal fallback.
%
% INPUTS (the benchmark planar four-cable setup)
%   WDM.param.M       2-by-4 anchor positions [mm], in fixed CCW cable order.
%   WDM.param.lim_inf  4-by-1 lower tension bounds [N]; all equal and positive.
%   WDM.param.lim_sup  4-by-1 upper tension bounds [N]; all equal, above lim_inf.
%   WDM.param.TOLL     Numerical zero/geometric cutoff (benchmark value 1e-14).
%   WDM.verbose       Logical scalar; use false for numerical calls.
%   P                 2-by-1 known platform position [mm].
%   F                 2-by-1 requested Cartesian force [N].
% W has unit directions from P to the anchors; the sign convention is W*T=F.
% Return a new output struct from a fresh parameter struct on each call.
% See docs/FUNCTION_REFERENCE.txt for the complete common contract and checks.
%
% Requires Optimization Toolbox. Geometry and H=I are rebuilt on each call.
% The constant options object is cached per MATLAB process; x0=[] and no
% previous tension solution is reused. Benchmark options: interior-point-convex,
% ConstraintTolerance=1e-10, OptimalityTolerance=1e-10, MaxIterations=200,
% Display=off. Other quadprog options use the installed MATLAB defaults.
%
% OUTPUTS (added to the returned WDM)
%   T                Solution/last quadprog iterate [N]; may be empty on failure.
%   objective        Scalar 0.5*T'*T [N^2], as returned by quadprog.
%   exitflag         Native quadprog status; positive is required for success.
%   iter             Native iteration count.
%   resid            Native firstorderopt when exposed; otherwise NaN.
%   primal_resid     norm(W*T-F,Inf) [N], if T is finite and correctly sized.
%   normT, errT      norm(T,2), norm(W*T-F,2) [N].
%   bound_violation  max([0;lim_inf-T;T-lim_sup]) [N].
%   case             NaN.
%   scale_dir        0 only for positive exitflag, bound violation<=1e-9 N
%                    and force error<=1e-8 N; otherwise -1.
% Diagnostic scalars remain NaN without a finite correctly sized T. Underlying
% solver errors may propagate. See docs/METHOD_SOURCES.txt for provenance.

persistent options
if isempty(options)
    options=optimoptions('quadprog','Algorithm','interior-point-convex', ...
        'Display','off','ConstraintTolerance',1e-10, ...
        'OptimalityTolerance',1e-10,'MaxIterations',200);
end
W=WDM.param.M-P;
for cable=1:size(W,2)
    length_cable=norm(W(:,cable));
    if length_cable>WDM.param.TOLL
        W(:,cable)=W(:,cable)/length_cable;
    else
        % Retain the benchmark's zero-length-cable convention. Such poses
        % are reported separately from the strictly interior domain.
        W(:,cable)=0;
    end
end
[T,objective,exitflag,details]=quadprog(eye(size(W,2)),[],[],[],W,F, ...
    WDM.param.lim_inf,WDM.param.lim_sup,[],options);
WDM.T=T;
WDM.objective=objective;
WDM.exitflag=exitflag;
WDM.iter=details.iterations;
WDM.case=NaN;
WDM.resid=NaN;
if isfield(details,'firstorderopt'),WDM.resid=details.firstorderopt;end
WDM.primal_resid=NaN;
WDM.normT=NaN;
WDM.errT=NaN;
WDM.bound_violation=NaN;
WDM.scale_dir=-1;
if numel(T)==size(W,2) && all(isfinite(T))
    WDM.normT=norm(T);
    WDM.errT=norm(W*T-F);
    WDM.primal_resid=norm(W*T-F,Inf);
    WDM.bound_violation=max([0;WDM.param.lim_inf-T;T-WDM.param.lim_sup]);
    if exitflag>0 && WDM.bound_violation<=1e-9 && WDM.errT<=1e-8
        WDM.scale_dir=0;
    end
end
end
