function WDM=wt_pott_v2(WDM,P,F)
% WT_POTT_V2 Improved closed-form cable-tension distribution (ICFM).
%   OUT = wt_pott_v2(WDM,P,F) starts at the midpoint tension vector, projects
%   with MATLAB pinv, and clamps/removes violating cables in up to three
%   nested steps. Lower-bound violations are tested before upper violations.
%   This is not a quadprog call and is not an L2-optimality certificate.
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
% The common bounds are read from the first vector entries during branching;
% use equal bounds across cables. No WEC or force adjustment occurs here.
%
% OUTPUTS (added to the returned WDM)
%   T          4-by-1 actual tensions [N], in original cable order.
%   normT      Scalar norm(T,2) [N].
%   errT       Scalar norm(W*T-F,2) [N].
%   scale_dir  1 if the final T violates a bound; otherwise 0.
%   active     4-by-1 logical test T>lim_inf.
%   cable_n    Number with abs((T-lim_inf)./lim_inf)>1e-6.
%   case       NaN (no Direct-method configuration).
%   param.M_eff Copy of the supplied anchor matrix.
% Check the result independently; scale_dir does not certify equilibrium.
% No custom pinv tolerance, iterative stopping rule, or warm-start state is
% used. Method reference: A. Pott (2014), manuscript reference [27]; see
% METHOD_SOURCES.txt. The captured file does not establish author-code origin.

Tmin=WDM.param.lim_inf;
Tmax=WDM.param.lim_sup;
WDM.scale_dir=0;
WDM.param.M_eff=WDM.param.M;
% WDM.str_align=[];

WDM.param.M_eff=WDM.param.M;

%calcolo versori
v = WDM.param.M_eff-repmat(P,1,size(WDM.param.M_eff,2));
for i=1:size(WDM.param.M_eff,2)
    d=norm(v(:,i)); d=(d>WDM.param.TOLL)*d;
    if d==0 
        v(:,i)=[0 0]';
    else
        v(:,i)=v(:,i)/d;
    end
end

v_pott = v;
F_pott = -F;
T_mean = (Tmin + Tmax) / 2;
% T_mean = Tmin;
T_pott = T_mean - pinv(v_pott) * (F_pott + v_pott * T_mean);
Tmin_pott = Tmin(1);
Tmax_pott = Tmax(1);


if any(T_pott < Tmin_pott)
    [T_pott, v_pott, T_mean, F_pott, index_1] = pott_less(T_pott, v_pott, F_pott, T_mean, Tmin_pott);
    
    if any(T_pott < Tmin_pott)
       [T_pott, v_pott, T_mean, F_pott, index_2] = pott_less(T_pott, v_pott, F_pott, T_mean, Tmin_pott);
       if any(T_pott < Tmin_pott)
           [T_pott, v_pott, T_mean, F_pott, index] = pott_less(T_pott, v_pott, F_pott, T_mean, Tmin_pott);
           T_pott = [T_pott(1:index-1); Tmin_pott; T_pott(index:end)];
        elseif any(T_pott > Tmax_pott)
            [T_pott, v_pott, T_mean, F_pott, index] = pott_greater(T_pott, v_pott, F_pott, T_mean, Tmax_pott);
            T_pott = [T_pott(1:index-1); Tmax_pott; T_pott(index:end)];
        end
       T_pott = [T_pott(1:index_2-1); Tmin_pott; T_pott(index_2:end)];
    elseif any(T_pott > Tmax_pott)
        [T_pott, v_pott, T_mean, F_pott, index_2] = pott_greater(T_pott, v_pott, F_pott, T_mean, Tmax_pott);
        if any(T_pott < Tmin_pott)
           [T_pott, v_pott, T_mean, F_pott, index] = pott_less(T_pott, v_pott, F_pott, T_mean, Tmin_pott);
           T_pott = [T_pott(1:index-1); Tmin_pott; T_pott(index:end)];
        elseif any(T_pott > Tmax_pott)
            [T_pott, v_pott, T_mean, F_pott, index] = pott_greater(T_pott, v_pott, F_pott, T_mean, Tmax_pott);
            T_pott = [T_pott(1:index_2-1); Tmax_pott; T_pott(index_2:end)];
        end
        T_pott = [T_pott(1:index_2-1); Tmax_pott; T_pott(index_2:end)];
    end

    T_pott = [T_pott(1:index_1-1); Tmin(index_1); T_pott(index_1:end)];

elseif any(T_pott > Tmax_pott)
    [T_pott, v_pott, T_mean, F_pott, index_1] = pott_greater(T_pott, v_pott, F_pott, T_mean, Tmax_pott);
    if any(T_pott < Tmin_pott)
        [T_pott, v_pott, T_mean, F_pott, index_2] = pott_less(T_pott, v_pott, F_pott, T_mean, Tmin_pott);
        if any(T_pott < Tmin_pott)
           [T_pott, v_pott, T_mean, F_pott, index] = pott_less(T_pott, v_pott, F_pott, T_mean, Tmin_pott);
           T_pott = [T_pott(1:index-1); Tmin_pott; T_pott(index:end)];
        elseif any(T_pott > Tmax_pott)
            [T_pott, v_pott, T_mean, F_pott, index] = pott_greater(T_pott, v_pott, F_pott, T_mean, Tmax_pott);
            T_pott = [T_pott(1:index-1); Tmax_pott; T_pott(index:end)];
        end
        T_pott = [T_pott(1:index_2-1); Tmin_pott; T_pott(index_2:end)];
    elseif any(T_pott > Tmax_pott)
        [T_pott, v_pott, T_mean, F_pott, index_2] = pott_greater(T_pott, v_pott, F_pott, T_mean, Tmax_pott);
        if any(T_pott < Tmin_pott)
           [T_pott, v_pott, T_mean, F_pott, index] = pott_less(T_pott, v_pott, F_pott, T_mean, Tmin_pott);
           T_pott = [T_pott(1:index-1); Tmin_pott; T_pott(index:end)];
        elseif any(T_pott > Tmax_pott)
            [T_pott, v_pott, T_mean, F_pott, index] = pott_greater(T_pott, v_pott, F_pott, T_mean, Tmax_pott);
            T_pott = [T_pott(1:index-1); Tmax_pott; T_pott(index:end)];
        end
        T_pott = [T_pott(1:index_2-1); Tmax_pott; T_pott(index_2:end)];
    end

    T_pott = [T_pott(1:index_1-1); Tmax(index_1); T_pott(index_1:end)];
end


T = T_pott;
%% QP
% % Forza=WDM.F';
% Forza=F';
% T=quadprog(eye(size(WDM.param.M_eff,2)),[],[],[],v,Forza,Tmin,[],Tmin,opt)
%% output
WDM.T=T;
if (sum(WDM.T<=Tmax)~=length(WDM.T) || sum(WDM.T>=Tmin)~=length(WDM.T))
    WDM.scale_dir=1;
end
WDM.normT=sqrt(WDM.T'*WDM.T);
WDM.errT=norm(v*WDM.T-F);
WDM.cable_n=sum(abs((WDM.T-Tmin)./Tmin)>1e-6);
WDM.active=WDM.T>Tmin;
WDM.case=nan;

% global T_pott_saved T_qp_saved
% T_pott_saved = [T_pott_saved, T_pott];
% T_qp_saved = [T_qp_saved, T];



end
function [T_pott, v_pott, T_mean, F_pott, index] = pott_less(T_pott, v_pott, F_pott, T_mean, Tmin_pott)
    [~, index] = min(T_pott);
    F_pott = Tmin_pott * v_pott(:,index) + F_pott;
    v_pott(:,index) = [];
    T_mean(index) = [];
    T_pott = T_mean - pinv(v_pott) * (F_pott + v_pott * T_mean);
end
function [T_pott, v_pott, T_mean, F_pott, index] = pott_greater(T_pott, v_pott, F_pott, T_mean, Tmax_pott)
    [~, index] = max(T_pott);
    F_pott = Tmax_pott * v_pott(:,index) + F_pott;
    v_pott(:,index) = [];
    T_mean(index) = [];
    T_pott = T_mean - pinv(v_pott) * (F_pott + v_pott * T_mean);
end
