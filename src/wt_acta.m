function WDM=wt_acta(WDM,P,F)
% WT_ACTA Analytic-centre cable-tension distribution.
%   OUT = wt_acta(WDM,P,F) minimizes the unit strict logarithmic barrier
%   -sum(log(T-lo)+log(hi-T)) subject to W*T=F. A finite analytic centre
%   requires a strictly feasible tension vector and full-row-rank W.
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
% Each call starts at T=(lo+hi)/2, with zero equality multipliers. WDM.acta
% may override tol, maxit, alpha, beta, maxls. Defaults: tol=1e-10 on the
% ABSOLUTE combined KKT 2-norm; maxit=50 accepted Newton updates; alpha=.01,
% beta=.5; maxls=60 shared step reductions for domain/decrease checks.
% The Newton system is reduced to a 2-by-2 Schur complement. No previous
% solution, WEC, or constant-options cache is used.
%
% OUTPUTS (added to the returned WDM)
%   T          4-by-1 last tension iterate [N], including on failure.
%   normT      norm(T,2) [N].
%   errT       norm(W*T-F,2) [N].
%   iter       Accepted Newton updates.
%   resid      Final absolute norm([grad(phi)+W'*lambda;W*T-F],2).
%   exitflag   1: residual<=tol; 0: budget exhausted; -1: singular-system
%              safeguard; -3: line-search safeguard.
%   scale_dir  0: converged and within bounds; -1: not converged; 1: bounds.
%   case       NaN.
% Flag 0 cannot distinguish an infeasible target from a difficult feasible
% cold start. Independently check the returned tensions and force residual.
%
% Equation-based MATLAB implementation; no original-author code was used.
% V. Di Paola, A. Goldsztejn, M. Zoppi and S. Caro (2024), Analytic Centre
% Based Tension Distribution for Cable-Driven Platforms (CDPs), Journal of
% Mechanisms and Robotics 16(8), 081018. DOI: 10.1115/1.4065244.
% Section 4, Eqs. 12-18 and Algorithm 1; infeasible-start Newton follows
% Boyd and Vandenberghe, Convex Optimization (2004), Algorithm 10.2.
% Defaults and cold-start policy are implementation choices; see docs.

Tmin = WDM.param.lim_inf;
Tmax = WDM.param.lim_sup;
TOL = 1e-10; MAXIT = 50; ALPHA = 0.01; BETA = 0.5; MAXLS = 60;
if isfield(WDM,'acta')
    o = WDM.acta;
    if isfield(o,'tol'),   TOL   = o.tol;   end
    if isfield(o,'maxit'), MAXIT = o.maxit; end
    if isfield(o,'alpha'), ALPHA = o.alpha; end
    if isfield(o,'beta'),  BETA  = o.beta;  end
    if isfield(o,'maxls'), MAXLS = o.maxls; end
end

% ---- unit cable vectors (platform -> motor), same as the other methods ----
v = WDM.param.M - repmat(P,1,size(WDM.param.M,2));
for i=1:size(WDM.param.M,2)
    d=norm(v(:,i)); d=(d>WDM.param.TOLL)*d;
    if d==0, v(:,i)=[0;0]; else, v(:,i)=v(:,i)/d; end
end
vt = v';

% ---- initial iterate: centre of the tension box, zero multipliers ----
tau = 0.5*(Tmin+Tmax);
mu  = zeros(size(v,1),1);
d1  = tau - Tmin;  d2 = Tmax - tau;
rd  = (1./d2 - 1./d1) + vt*mu;      % grad phi + W'*mu
rp  = v*tau - F;                    % W*tau + we
nr  = sqrt(rd'*rd + rp'*rp);        % ||Eq.(15)||_2

iter = 0; exitflag = 0;
while nr > TOL && iter < MAXIT
    % ---- Newton step (Eq. 16) by block elimination ----
    hi = 1./(1./(d1.*d1) + 1./(d2.*d2));      % H^-1 (diagonal of the Hessian inverse)
    Wh = v.*hi';                              % W*H^-1
    S  = Wh*vt;                               % W*H^-1*W'  (2x2)
    b  = rp - Wh*rd;
    dS = S(1,1)*S(2,2) - S(1,2)*S(2,1);
    if ~(dS > 1e-12*S(1,1)*S(2,2))            % W rank deficient -> Eq.(16) singular
        exitflag = -1; break
    end
    dmu  = [S(2,2)*b(1) - S(1,2)*b(2); S(1,1)*b(2) - S(2,1)*b(1)]/dS;
    dtau = -hi.*(rd + vt*dmu);

    % ---- backtracking line search (B&V Alg. 10.2) ----
    t = 1; nls = 0;
    taun = tau + dtau; d1n = taun - Tmin; d2n = Tmax - taun;
    while any(d1n <= 0) || any(d2n <= 0)      % keep tau strictly inside the box
        nls = nls + 1;
        if nls > MAXLS, break; end
        t = BETA*t;
        taun = tau + t*dtau; d1n = taun - Tmin; d2n = Tmax - taun;
    end
    if nls > MAXLS, exitflag = -3; break; end
    mun = mu + t*dmu;
    rdn = (1./d2n - 1./d1n) + vt*mun;
    rpn = v*taun - F;
    nrn = sqrt(rdn'*rdn + rpn'*rpn);
    while ~(nrn <= (1-ALPHA*t)*nr)            % sufficient decrease of ||r||_2
        nls = nls + 1;
        if nls > MAXLS, break; end
        t = BETA*t;
        taun = tau + t*dtau; d1n = taun - Tmin; d2n = Tmax - taun;
        mun = mu + t*dmu;
        rdn = (1./d2n - 1./d1n) + vt*mun;
        rpn = v*taun - F;
        nrn = sqrt(rdn'*rdn + rpn'*rpn);
    end
    if nls > MAXLS, exitflag = -3; break; end

    % ---- update ----
    tau = taun; mu = mun; d1 = d1n; d2 = d2n; rd = rdn; rp = rpn; nr = nrn;
    iter = iter + 1;
end
if nr <= TOL, exitflag = 1; end

% ---- outputs ----
WDM.T        = tau;
WDM.normT    = sqrt(tau'*tau);
WDM.errT     = norm(v*tau - F);
WDM.iter     = iter;
WDM.resid    = nr;
WDM.exitflag = exitflag;
if exitflag < 1
    WDM.scale_dir = -1;
elseif sum(tau<=Tmax)~=length(tau) || sum(tau>=Tmin)~=length(tau)
    WDM.scale_dir = 1;
else
    WDM.scale_dir = 0;
end
WDM.case = nan;
end
