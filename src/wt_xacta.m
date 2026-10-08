function WDM = wt_xacta(WDM,P,F)
% WT_XACTA Extended analytic-centre cable-tension distribution.
%   OUT = wt_xacta(WDM,P,F) runs relaxed analytic centring (RAC), then EAC
%   when activation requires it. It accepts the raw requested force and
%   does not call WEC. EAC permits force error while keeping strict bounds.
%
% INPUTS (the reviewed planar four-cable setup)
%   WDM.param.M       2-by-4 anchor positions [mm], in fixed CCW cable order.
%   WDM.param.lim_inf  4-by-1 lower tension bounds [N]; all equal and positive.
%   WDM.param.lim_sup  4-by-1 upper tension bounds [N]; all equal, above lim_inf.
%   WDM.param.TOLL     Numerical zero/geometric cutoff (reviewed value 1e-14).
%   WDM.verbose       Logical scalar; use false for numerical calls.
%   P                 2-by-1 known platform position [mm].
%   F                 2-by-1 requested Cartesian force [N].
% W has unit directions from P to the anchors; the sign convention is W*T=F.
% Return a new output struct from a fresh parameter struct on each call.
% See docs/FUNCTION_REFERENCE.txt for the complete common contract and checks.
%
% The reviewed contract is planar 2-by-4; this function also validates general
% n-by-m anchors (n<=m), n-vector P/F, and scalar or m-vector ordered bounds.
% RAC starts from the affine projection of the box midpoint, with zero null
% coordinates; its relaxed barrier allows intermediate bound violations.
% EAC reuses this call's RAC iterate if interior, otherwise its shrunken-box
% projection. No final clipping or previous-solution warm start is performed.
% A one-entry exact-key cache retains validated constants/options only; each
% call validates P/F, rebuilds geometry and creates its initial iterate.
%
% WDM.xacta options and defaults:
%   delta=.1 N, gamma=500, eta=1e-5, k2=k3=1, tol=1e-9, force_tol=1e-8 N;
%   maxit=100 per stage; rac_maxit/eac_maxit optionally override it;
%   alpha=.01, beta=.5, maxls=60 (61 trial indices, 0:60).
% Both stages require relative equality residual<=tol and absolute infinity
% norm<=force_tol. Relative stationarity uses
%   max(tol,min(sqrt(eps),4*roundoff_floor)).
% The returned floor/tolerance make the precision allowance explicit. In EAC,
% the equality residual is augmented: W*T-F-c*lambda, not raw force error.
% See FUNCTION_REFERENCE for the exact scales and activation/slack equations.
%
% OUTPUTS (added to the returned WDM)
%   T, T_rac          Final/RAC m-by-1 actual tensions [N].
%   normT, errT       norm(T,2), norm(W*T-F,2) [N].
%   achieved_force    n-by-1 W*T [N].
%   mu, lambda, s     Activation scalar, n multipliers, n force slacks.
%   stage            1: RAC final; 2: EAC final.
%   exitflag         1: convergence; 0: budget; -1: rank/factorization;
%                    -2: nonfinite state; -3: line search cannot progress.
%   iter             Sum of accepted updates in both stages.
%   scale_dir, case  0 on success, -1 otherwise; case=NaN.
% Residuals, stage flags/iterations and work counters are listed in the function
% reference. Unrun EAC residual/flag fields are NaN; counters are zero. Invalid
% inputs/options raise errors. Native convergence does not imply W*T=F in EAC.
%
% Equation-based MATLAB implementation; no original-author code was used.
% D. Dona', V. Di Paola, A. Trevisani and M. Zoppi (2026), X-ACTA: eXtended
% Analytic Center Tension distribution Algorithm for fixed and mobile
% cable-driven-parallel-robot, arXiv:2607.08265v1, Eqs. 14-21.
% Paper source recorded by the measured header:
% https://arxiv.org/html/2607.08265v1
% The numerical formulations and safeguards are implementation choices.

[W,lo,hi,F,o] = inputs(WDM,P,F);
n = size(W,1); m = size(W,2);
tau0 = (lo+hi)/2;
lambda0 = zeros(n,1);

% The full-rank assumption is checked before either optimization stage.
WW = W*W';
if rcond(WW) < 1e-14
    rac = empty_state(tau0,lambda0,-1);
else
    rac = rac_stage(W,WW,F,lo,hi,tau0,o);
end
tau = rac.tau; lambda = rac.lambda;
mu = NaN; s = zeros(n,1); stage = 1;
eac = empty_state(nan(m,1),nan(n,1),NaN);
if rac.flag == 1
    shrunk = min(max(rac.tau,lo+o.delta),hi-o.delta);
    distance = abs(rac.tau-shrunk);
    mu = tanh(sum(distance.^(o.k3+1)));
    if any(distance > 0)
        stage = 2;
        setup_evals = 0;
        % Reuse this call's RAC solution where the true log is defined.
        if all(rac.tau > lo & rac.tau < hi)
            tau0 = rac.tau;
            lambda0 = rac.lambda;
        else
            tau0 = shrunk;
            [g,~] = barrier(tau0,lo,hi,false,o.delta,o.k2);
            lambda0 = -WW\(W*g);
            setup_evals = 1;
        end
        a = o.eta+o.gamma*mu^2;
        c = mu^2/(2*a);
        eac = newton_stage(W,F,lo,hi,tau0,lambda0,c,false,o,o.eac_maxit);
        eac.barrier_evals = eac.barrier_evals+setup_evals;
        tau = eac.tau; lambda = eac.lambda;
        s = -mu*lambda/(2*a);
    end
end

if stage == 1, final = rac; else, final = eac; end
WDM.T = tau;
WDM.normT = norm(tau);
WDM.achieved_force = W*tau;
WDM.errT = norm(WDM.achieved_force-F);
WDM.T_rac = rac.tau;
WDM.rac_lambda = rac.lambda;
WDM.mu = mu;
WDM.lambda = lambda;
WDM.s = s;
WDM.stage = stage;
WDM.iter = rac.iter+eac.iter;
WDM.rac_iters = rac.iter;
WDM.eac_iters = eac.iter;
WDM.rac_iter = rac.iter;
WDM.eac_iter = eac.iter;
WDM.exitflag = final.flag;
WDM.rac_exitflag = rac.flag;
WDM.eac_exitflag = eac.flag;
WDM.resid = final.scaled;
WDM.raw_resid = final.raw;
WDM.primal_resid = final.primal;
WDM.stationarity_resid = final.stationarity;
WDM.stationarity_relative = final.dual_relative;
WDM.stationarity_tolerance = final.stationarity_tolerance;
WDM.roundoff_floor = final.roundoff_floor;
WDM.rac_resid = rac.scaled;
WDM.eac_resid = eac.scaled;
WDM.rac_barrier_evals = rac.barrier_evals;
WDM.eac_barrier_evals = eac.barrier_evals;
WDM.barrier_evals = rac.barrier_evals+eac.barrier_evals;
WDM.rac_backtracks = rac.backtracks;
WDM.eac_backtracks = eac.backtracks;
WDM.backtracks = rac.backtracks+eac.backtracks;
WDM.scale_dir = 0;
if final.flag ~= 1, WDM.scale_dir = -1; end
WDM.case = NaN;
end

function out = rac_stage(W,WW,F,lo,hi,tau0,o)
% Feasible affine coordinates make this ordinary unconstrained Newton.
base = tau0+W'*(WW\(F-W*tau0));
Z = null(W);
z = zeros(size(Z,2),1);
tau = base;
iter = 0; flag = 0;
force_scale = max(1,norm(F,Inf));
[g,h,f] = barrier(tau,lo,hi,true,o.delta,o.k2);
barrier_evals = 1; backtracks = 0;
while true
    lambda = -(WW\(W*g));
    Wlambda = W'*lambda;
    rd = g+Wlambda; rp = W*tau-F;
    primal = norm(rp,Inf);
    dual_scale = max([1;abs(g);abs(Wlambda)]);
    dual_relative = norm(rd,Inf)/dual_scale;
    roundoff_floor = norm(h.*eps(abs(tau)),Inf)/dual_scale;
    stationarity_tolerance = max(o.tol,min(sqrt(eps),4*roundoff_floor));
    scaled = max(dual_relative,primal/force_scale);
    if ~all(isfinite([tau;g;h;lambda;scaled])) || any(h <= 0)
        flag = -2; break
    end
    if dual_relative <= stationarity_tolerance && primal/force_scale <= o.tol && primal <= o.force_tol
        flag = 1; break
    end
    if iter >= o.rac_maxit, break; end
    Q = Z'*(h.*Z);
    [L,p] = chol((Q+Q')/2,'lower');
    if p ~= 0, flag = -1; break; end
    dz = -(L'\(L\(Z'*g)));
    dtau = Z*dz;
    slope = g'*dtau;
    accepted = false; step = 1;
    for ls = 0:o.maxls
        zn = z+step*dz;
        taun = base+Z*zn;
        [gn,hn,fn] = barrier(taun,lo,hi,true,o.delta,o.k2);
        barrier_evals = barrier_evals+1;
        if scaled <= 1e-4
            lambdan = -(WW\(W*gn));
            rdn = gn+W'*lambdan;
            acceptable = norm(rdn) <= (1-o.alpha*step)*norm(rd);
        else
            acceptable = fn <= f+o.alpha*step*slope+8*eps(max(1,abs(f)));
        end
        % Roundoff allowance affects line search only; success still
        % requires the independently recomputed KKT and primal residuals.
        if isfinite(fn) && acceptable
            accepted = true; break
        end
        backtracks = backtracks+1;
        step = o.beta*step;
    end
    if ~accepted || all(taun == tau), flag = -3; break; end
    % These are already the derivatives at the accepted point. Reusing
    % them preserves the iterates and avoids a second barrier evaluation.
    tau = taun; z = zn; g = gn; h = hn; f = fn; iter = iter+1;
end
out.tau = tau; out.lambda = lambda; out.flag = flag; out.iter = iter;
out.scaled = scaled; out.raw = norm([rd;rp]);
out.primal = primal; out.stationarity = norm(rd,Inf);
out.dual_relative = dual_relative;
out.stationarity_tolerance = stationarity_tolerance;
out.roundoff_floor = roundoff_floor;
out.barrier_evals = barrier_evals; out.backtracks = backtracks;
end

function out = newton_stage(W,F,lo,hi,tau,lambda,c,relaxed,o,maxit)
% Residual merit uses FIXED scales, preserving its Newton derivative -r.
[g,h] = barrier(tau,lo,hi,relaxed,o.delta,o.k2);
Wlambda = W'*lambda;
dual_scale = max([1; abs(g); abs(Wlambda)]);
force_scale = max(1,norm(F,Inf));
rd = g+Wlambda;
rp = W*tau-F-c*lambda;
merit = norm([rd/dual_scale;rp/force_scale]);
iter = 0; flag = 0;
barrier_evals = 1; backtracks = 0;
cI = c*eye(size(W,1));
while true
    current_dual_scale = max([1;abs(g);abs(Wlambda)]);
    dual_relative = norm(rd,Inf)/current_dual_scale;
    roundoff_floor = norm(h.*eps(abs(tau)),Inf)/current_dual_scale;
    stationarity_tolerance = max(o.tol,min(sqrt(eps),4*roundoff_floor));
    primal = norm(rp,Inf);
    scaled = max(dual_relative,primal/force_scale);
    if ~all(isfinite([tau;lambda;g;h;rd;rp;scaled])) || any(h <= 0)
        flag = -2; break
    end
    if dual_relative <= stationarity_tolerance && primal/force_scale <= o.tol && primal <= o.force_tol
        flag = 1; break
    end
    if iter >= maxit, break; end
    hinv = 1./h;
    Wh = W.*hinv';
    S = Wh*W'+cI;
    % Symmetrize the small Schur system against accumulation roundoff.
    S = (S+S')/2;
    [L,p] = chol(S,'lower');
    if p ~= 0 || any(~isfinite(S(:)))
        flag = -1; break
    end
    % Solve directly for the multiplier at the full Newton step. This
    % equivalent form avoids cancellation of old and updated multipliers.
    lambda_target = L'\(L\(W*tau-F-Wh*g));
    dlambda = lambda_target-lambda;
    dtau = -hinv.*(g+W'*lambda_target);
    if any(~isfinite([dtau;dlambda]))
        flag = -2; break
    end

    step = 1;
    if ~relaxed
        % Fraction-to-boundary is a numerical domain safeguard, not a
        % clipping operation on the optimizer's returned tensions.
        neg = dtau < 0; pos = dtau > 0;
        if any(neg), step = min(step,.99*min((lo(neg)-tau(neg))./dtau(neg))); end
        if any(pos), step = min(step,.99*min((hi(pos)-tau(pos))./dtau(pos))); end
    end
    accepted = false;
    e = W*tau-F;
    objective_slope = c*(g'*dtau)+e'*(W*dtau);
    for ls = 0:o.maxls
        taun = tau+step*dtau;
        lambdan = lambda+step*dlambda;
        if relaxed || all(taun > lo & taun < hi)
            [gn,hn] = barrier(taun,lo,hi,relaxed,o.delta,o.k2);
            barrier_evals = barrier_evals+1;
            Wlambdan = W'*lambdan;
            rdn = gn+Wlambdan;
            rpn = W*taun-F-c*lambdan;
            meritn = norm([rdn/dual_scale;rpn/force_scale]);
            if c > 1e-10 && scaled > 1e-4
                % Accurate log and squared-error objective differences;
                % there is no eta/mu^2 division or enormous penalty weight.
                change_barrier = -sum(log1p((taun-tau)./(tau-lo))) ...
                    -sum(log1p(-(taun-tau)./(hi-tau)));
                en = W*taun-F;
                change = c*change_barrier+.5*(en-e)'*(en+e);
                allowance = 8*eps(max(realmin,.5*(e'*e)));
                acceptable = change <= o.alpha*step*objective_slope+allowance;
            else
                acceptable = meritn <= (1-o.alpha*step)*merit;
            end
            if isfinite(meritn) && acceptable
                accepted = true; break
            end
        end
        backtracks = backtracks+1;
        step = o.beta*step;
    end
    if ~accepted || all(taun == tau) && all(lambdan == lambda)
        flag = -3; break
    end
    tau = taun; lambda = lambdan; g = gn; h = hn;
    Wlambda = Wlambdan;
    rd = rdn; rp = rpn; merit = meritn;
    iter = iter+1;
end
out.tau = tau; out.lambda = lambda; out.flag = flag; out.iter = iter;
out.scaled = scaled;
out.raw = norm([rd;rp]);
out.primal = norm(rp,Inf);
out.stationarity = norm(rd,Inf);
out.dual_relative = dual_relative;
out.stationarity_tolerance = stationarity_tolerance;
out.roundoff_floor = roundoff_floor;
out.barrier_evals = barrier_evals; out.backtracks = backtracks;
end

function [g,h,f] = barrier(tau,lo,hi,relaxed,delta,k2)
d1 = tau-lo; d2 = hi-tau;
if relaxed
    [g1,h1,f1] = relaxed_log(d1,delta,k2);
    [g2,h2,f2] = relaxed_log(d2,delta,k2);
    g = g1-g2; h = h1+h2;
    f = sum(f1+f2);
else
    g = -1./d1+1./d2;
    h = 1./d1.^2+1./d2.^2;
    % EAC uses accurate log1p differences in its line search, so it needs
    % derivatives here but never the absolute logarithmic objective.
    if nargout > 2, f = -sum(log(d1)+log(d2)); end
end
end

function [g,h,f] = relaxed_log(d,delta,k2)
inside = d >= delta;
g = zeros(size(d)); h = g; f = g;
g(inside) = -1./d(inside);
h(inside) = 1./d(inside).^2;
f(inside) = -log(d(inside));
outside = ~inside;
if any(outside)
    z = (delta-d(outside))/delta;
    if k2 == 1
        g(outside) = -(1+z)/delta;
        h(outside) = 1/delta^2;
        f(outside) = -log(delta)+z+.5*z.^2;
    else
        % Eq. 15 expressed with z>=0: all polynomial terms have positive
        % coefficients; its negative first and positive second derivatives
        % avoid alternating sums and retain convexity for every k2>=1.
        gp = ones(size(z)); hp = zeros(size(z)); fp = -log(delta)+z;
        for j = 2:k2+1
            gp = gp+z.^(j-1);
            hp = hp+(j-1)*z.^(j-2);
            fp = fp+z.^j/j;
        end
        g(outside) = -gp/delta;
        h(outside) = hp/delta^2;
        f(outside) = fp;
    end
end
end

function out = empty_state(tau,lambda,flag)
out.tau = tau; out.lambda = lambda; out.flag = flag; out.iter = 0;
out.scaled = NaN; out.raw = NaN; out.primal = NaN; out.stationarity = NaN;
out.dual_relative = NaN; out.stationarity_tolerance = NaN; out.roundoff_floor = NaN;
out.barrier_evals = 0; out.backtracks = 0;
end

function [W,lo,hi,F,o] = inputs(WDM,P,F)
% Cache only successfully checked constants; no pose, force, iterate or
% factorization is retained. Each MATLAB worker owns its own one-entry cache.
persistent cache
if ~isstruct(WDM) || ~isscalar(WDM) || ~isfield(WDM,'param') ...
        || ~isstruct(WDM.param) || ~isscalar(WDM.param)
    error('wt_xacta:InvalidInput','WDM.param must be a scalar parameter struct.');
end
p = WDM.param;
required = {'M','TOLL','lim_inf','lim_sup'};
for j = 1:numel(required)
    if ~isfield(p,required{j})
        error('wt_xacta:InvalidInput','Missing WDM.param.%s.',required{j});
    end
end
has_options = isfield(WDM,'xacta');
raw_options = [];
if has_options
    raw_options = WDM.xacta;
    if ~isstruct(raw_options) || ~isscalar(raw_options)
        error('wt_xacta:InvalidOption','WDM.xacta must be a scalar options struct.');
    end
end
if ~isempty(cache) && same_configuration(p,raw_options,has_options,cache)
    lo = cache.lo; hi = cache.hi; o = cache.o; n = cache.n;
else
    [lo,hi,o,n] = validate_configuration(p,raw_options,has_options);
    % Assign the cache only after the entire constant configuration passes.
    fresh.M = p.M; fresh.TOLL = p.TOLL;
    fresh.lim_inf = p.lim_inf; fresh.lim_sup = p.lim_sup;
    fresh.raw_options = raw_options; fresh.has_options = has_options;
    fresh.lo = lo; fresh.hi = hi; fresh.o = o; fresh.n = n;
    if has_options
        fresh.option_names = fieldnames(raw_options);
    else
        fresh.option_names = {};
    end
    cache = fresh;
end
% Changing inputs are always validated, including on a configuration hit.
validateattributes(P,{'numeric'},{'real','finite','vector','numel',n},mfilename,'P');
validateattributes(F,{'numeric'},{'real','finite','vector','numel',n},mfilename,'F');
P = double(P(:)); F = double(F(:));
W = double(p.M)-P;
lengths = sqrt(sum(W.^2,1));
nonzero = lengths > p.TOLL;
W(:,nonzero) = W(:,nonzero)./lengths(nonzero);
W(:,~nonzero) = 0;
end

function same = same_configuration(p,raw_options,has_options,c)
% isequal compares every element and shape; explicit class comparisons also
% distinguish numerically equal integer/single/double configurations. No
% checksum or rounded key can cause stale options or bounds to be reused.
same = has_options == c.has_options ...
    && same_numeric(p.M,c.M) && same_numeric(p.TOLL,c.TOLL) ...
    && same_numeric(p.lim_inf,c.lim_inf) && same_numeric(p.lim_sup,c.lim_sup) ...
    && isequal(raw_options,c.raw_options);
if same && has_options
    for j = 1:numel(c.option_names)
        name = c.option_names{j};
        if ~same_numeric(raw_options.(name),c.raw_options.(name))
            same = false; return
        end
    end
end
end

function same = same_numeric(a,b)
% In particular, a complex array with zero imaginary parts must not bypass
% the real-input validation just because isequal regards its values equal.
same = strcmp(class(a),class(b)) && isreal(a) == isreal(b) ...
    && issparse(a) == issparse(b) && isequal(a,b);
end

function [lo,hi,o,n] = validate_configuration(p,raw_options,has_options)
% Full public-input validation runs for every new constant configuration.
validateattributes(p.M,{'numeric'},{'real','finite','2d','nonempty'},mfilename,'M');
[n,m] = size(p.M);
if n > m, error('wt_xacta:InvalidInput','M must have at least as many cables as force dimensions.'); end
validateattributes(p.TOLL,{'numeric'},{'real','finite','scalar','nonnegative'},mfilename,'TOLL');
validateattributes(p.lim_inf,{'numeric'},{'real','finite','nonempty'},mfilename,'lim_inf');
validateattributes(p.lim_sup,{'numeric'},{'real','finite','nonempty'},mfilename,'lim_sup');
lo = double(p.lim_inf(:)); hi = double(p.lim_sup(:));
if isscalar(lo), lo = repmat(lo,m,1); end
if isscalar(hi), hi = repmat(hi,m,1); end
if numel(lo) ~= m || numel(hi) ~= m || any(lo >= hi)
    error('wt_xacta:InvalidInput','Tension bounds must be ordered scalar values or one value per cable.');
end
o = struct('delta',.1,'gamma',500,'eta',1e-5,'k2',1,'k3',1, ...
    'tol',1e-9,'force_tol',1e-8,'maxit',100,'alpha',.01,'beta',.5,'maxls',60);
if has_options
    names = fieldnames(raw_options);
    supported = [fieldnames(o);{'rac_maxit';'eac_maxit'}];
    for j = 1:numel(names)
        if ~ismember(names{j},supported)
            error('wt_xacta:InvalidOption','Unknown X-ACTA option: %s.',names{j});
        end
        o.(names{j}) = raw_options.(names{j});
    end
end
if ~isfield(o,'rac_maxit'), o.rac_maxit = o.maxit; end
if ~isfield(o,'eac_maxit'), o.eac_maxit = o.maxit; end
positive = {'delta','gamma','eta','tol','force_tol'};
for j = 1:numel(positive)
    validateattributes(o.(positive{j}),{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,positive{j});
end
integers = {'k2','k3','maxit','rac_maxit','eac_maxit','maxls'};
for j = 1:numel(integers)
    validateattributes(o.(integers{j}),{'numeric'}, ...
        {'real','finite','scalar','integer','positive'},mfilename,integers{j});
end
validateattributes(o.alpha,{'numeric'},{'real','finite','scalar','>',0,'<',.5},mfilename,'alpha');
validateattributes(o.beta,{'numeric'},{'real','finite','scalar','>',0,'<',1},mfilename,'beta');
if 2*o.delta >= min(hi-lo)
    error('wt_xacta:InvalidOption','delta must be less than half of every tension range.');
end
end
