function WDM=wt_direct(WDM,P,F)
% WT_DIRECT Direct L2-optimal tension distribution.
%   OUT = wt_direct(WDM,P,F) computes a candidate with minimum
%   actual-tension norm for the four-cable geometry and feasible target of
%   the accompanying manuscript. WEC is a separate function; this routine
%   does not change F or call an iterative optimizer.
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
% Use the fixed CCW anchor order and a position strictly inside the convex
% anchor frame. Cyclic relabeling is not a general sort for unordered anchors
% or exterior/degenerate poses. The routine reads the first lower/
% upper bound as a common scalar; unequal cable bounds are unsupported here.
%
% OUTPUTS (added to the returned WDM)
%   T          4-by-1 actual tensions [N], in ORIGINAL cable order.
%   normT      Scalar norm(T,2) [N].
%   errT       Scalar norm(W*T-F,2) [N].
%   case       0: numerical r=0; 1..4: configurations I..IV; -1: saturation.
%   scale_dir  1: double-saturation/empty-interval flag; 0: no such flag.
% A zero flag is not a physical-feasibility certificate. A flagged branch
% retains its algebraic candidate; there is no fallback optimizer. Invalid
% or degenerate inputs may error or produce a nonfinite result. With verbose
% true, extra angle fields and N1dotN2 are computed by an optional 2-by-2 solve.
%
% METHOD
% Algorithm 1 of the accompanying manuscript: r=F/Tmin-W*ones(4,1), analytic
% tauP and N (Eqs. 6-7), full Pm=N'*N and q=N'*(tauP+1), axis intervals
% (Eqs. 11-12), scalar objective (Eq. 13), and clipping (Eq. 14).
% For nonzero r, alignment uses r/norm(r) directly, with four atan2 values,
% three ordering comparisons and seven sine evaluations. Configurations
% II/III use the alpha2/alpha1 axis respectively; IV compares both candidates.
% The zero-r branch returns the lower-bound tensions. The code has no
% previous-solution warm start. See docs/METHOD_SOURCES.txt for attribution.

Tmin    = WDM.param.lim_inf(1);
Tmax    = WDM.param.lim_sup(1);
Tau_max = (Tmax-Tmin)/Tmin;
TOL     = WDM.param.TOLL;

WDM.scale_dir = 0;
cable = [1;2;3;4]; N = zeros(4,2); TauP = zeros(4,1);

% ---- unit cable vectors (platform -> motor) ----
v = WDM.param.M - repmat(P,1,size(WDM.param.M,2));
for i=1:size(WDM.param.M,2)
    d=norm(v(:,i)); d=(d>TOL)*d;
    if d==0, v(:,i)=[0;0]; else, v(:,i)=v(:,i)/d; end
end

% Algorithm 1, line 2: normalize the requested force; numerical zero branch.
f  = F/Tmin;  r = f - sum(v,2);
nr = norm(r,2); nr = (nr>TOL)*nr;

if nr==0
    % r = 0  ->  tau = 0   (Case 0)
    Tau = zeros(4,1); casenum = 0;
else
    % Algorithm 1, lines 3-4: direct alignment and cyclic CCW relabeling.
    ur  = r/nr;
    Rrb = [ur(1) -ur(2); ur(2) ur(1)];
    vr  = Rrb'*v;
    Th  = atan2(vr(2,:),vr(1,:))';
    Th  = Th.*(Th>=0) + (Th+2*pi).*(Th<0);     % map to [0,2*pi)
    imax=1;
    if Th(imax)<Th(2), imax=2; end
    if Th(imax)<Th(3), imax=3; end
    if Th(imax)<Th(4), imax=4; end
    cable = [imax; sig(imax+1); sig(imax+2); sig(imax+3)];

    % Algorithm 1, lines 5-6; Eqs. (6)-(7): analytic tauP and N.
    s21 = sin(Th(cable(2))-Th(cable(1)));
    TauP = zeros(4,1);
    TauP(cable(1)) =  nr/s21*sin(Th(cable(2)));
    TauP(cable(2)) = -nr/s21*sin(Th(cable(1)));
    N = zeros(4,2);
    N(cable,:) = [ sin(Th(cable(3))-Th(cable(2)))  sin(Th(cable(4))-Th(cable(2)));
                  -sin(Th(cable(3))-Th(cable(1))) -sin(Th(cable(4))-Th(cable(1)));
                   s21                              0;
                   0                                s21 ];

    % Eq. (13): form full Pm=N' * N and q; only scalar minima are solved.
    q  = N'*(ones(4,1)+TauP);   % q = N'(tauP + 1)
    Pm = N'*N;                  % P = N'N

    alpha = [0;0];

    % ---- Corollary 1: double saturation -> unfeasible ----
    if (TauP(cable(1))>Tau_max) && (TauP(cable(2))>Tau_max)
        WDM.scale_dir = true; casenum = -1;
    else
        % configuration from signs:  N(cable(2),1) = -s31,  N(cable(1),2) = s42
        s31_le0 = N(cable(2),1) >= 0;     % -s31 >= 0  <=>  s31 <= 0
        s42_ge0 = N(cable(1),2) >= 0;     %  s42 >= 0

        if s31_le0 && s42_ge0
            % Configuration I : both columns of N >= 0, tauP dominates
            casenum = 1; alpha = [0;0];

        elseif s31_le0 && ~s42_ge0
            % Configuration II : optimum on the alpha2 axis
            casenum = 2;
            I2  = axis_interval(TauP, N(:,2), Tau_max);   % alpha2 axis (alpha1 = 0)
            a2u = -q(2)/Pm(2,2);
            if interval_empty(I2)
                WDM.scale_dir = true;
            else
                alpha = [0; clip(a2u,I2)];
            end

        elseif ~s31_le0 && s42_ge0
            % Configuration III : optimum on the alpha1 axis
            casenum = 3;
            I1  = axis_interval(TauP, N(:,1), Tau_max);   % alpha1 axis (alpha2 = 0)
            a1u = -q(1)/Pm(1,1);
            if interval_empty(I1)
                WDM.scale_dir = true;
            else
                alpha = [clip(a1u,I1); 0];
            end

        else
            % Configuration IV : take the better of the two axis candidates
            casenum = 4;
            I1  = axis_interval(TauP, N(:,1), Tau_max);   % alpha1 axis (alpha2 = 0)
            I2  = axis_interval(TauP, N(:,2), Tau_max);   % alpha2 axis (alpha1 = 0)
            a1u = -q(1)/Pm(1,1);
            a2u = -q(2)/Pm(2,2);
            haveI1 = ~interval_empty(I1);
            haveI2 = ~interval_empty(I2);
            if ~haveI1 && ~haveI2
                WDM.scale_dir = true;
            elseif haveI1 && ~haveI2
                alpha = [clip(a1u,I1); 0];
            elseif ~haveI1 && haveI2
                alpha = [0; clip(a2u,I2)];
            else
                aA = [clip(a1u,I1); 0];
                aB = [0; clip(a2u,I2)];
                fA = 0.5*aA'*Pm*aA + q'*aA;
                fB = 0.5*aB'*Pm*aB + q'*aB;
                if fA <= fB, alpha = aA; else, alpha = aB; end
            end
        end
    end

    Tau = TauP + N*alpha;
end

% Algorithm 1, line 22: T is in original cable order after indexed assignments.
WDM.T     = (ones(4,1)+Tau)*Tmin;
WDM.normT = sqrt(WDM.T'*WDM.T);
WDM.errT  = norm(v*WDM.T - F);
WDM.case  = casenum;

% analysis-only outputs (not part of Algorithm 1); computed only in verbose
% mode so the timed path stays free of the 2x2 solve and extra field writes
if WDM.verbose && casenum~=0
    qq = N(cable,:)'*(ones(4,1)+TauP(cable));
    alpha_opt_th = -(N(cable,:)'*N(cable,:))\qq;
    WDM.theta_q        = atan2(qq(2),qq(1));
    WDM.theta_alphamin = atan2(alpha_opt_th(2),alpha_opt_th(1));
    WDM.N1dotN2        = N(cable,1)'*N(cable,2);
end
end

%##########################################################################
% Local helper functions
%##########################################################################
function out=sig(t)        % cyclic cable index in 1..4
out = mod(t-1,4)+1;
end

function I=axis_interval(TauP,n,Tau_max)
% feasible interval  {t >= 0 : 0 <= TauP + t*n <= Tau_max}  returned as [lo hi]
lo=0; hi=Inf;
for i=1:4
    ni=n(i); tpi=TauP(i);
    if ni>0
        lo=max(lo,(0-tpi)/ni);
        hi=min(hi,(Tau_max-tpi)/ni);
    elseif ni<0
        hi=min(hi,(0-tpi)/ni);
        lo=max(lo,(Tau_max-tpi)/ni);
    else
        if tpi<-1e-12 || tpi>Tau_max+1e-12
            lo=Inf; hi=-Inf;     % constant tension out of bounds -> empty
        end
    end
end
I=[lo hi];
end

function e=interval_empty(I)
e = I(1) > I(2) + 1e-12;
end

function t=clip(x,I)
t = min(max(x,I(1)),I(2));
end
