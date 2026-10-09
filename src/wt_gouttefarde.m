function WDM=wt_gouttefarde(WDM,P,F)
% WT_GOUTTEFARDE Vertex-based cable-tension distribution (VTDA-L2).
%   OUT = wt_gouttefarde(WDM,P,F) constructs numerical null coordinates,
%   traces the feasible polygon, and tests vertex/edge L2 candidates. It may
%   fall back to the minimum-norm feasible vertex that it traced.
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
% This four-cable implementation uses null(W), an explicit 2-by-2 generalized
% inverse, and small linear solves. It has no rank regularization or WEC.
%
% OUTPUTS (when a new solution is produced)
%   T          4-by-1 actual tensions [N], in original cable order.
%   normT      Scalar norm(T,2) [N].
%   errT       Scalar norm(W*T-F,2) [N].
%   fallback   1: feasible traced-vertex fallback; 0: no such fallback.
%   scale_dir  1: detected bound violation (using TOLL); otherwise 0.
%   vertices   2-by-k traced vertices in numerical null coordinates.
%   param.M_eff Copy of the supplied anchor matrix.
% IMPORTANT: the walk stops after 100 processed visits and can return WITHOUT
% creating T; numerical errors may also interrupt the call. An incoming T is
% not cleared. Use a fresh parameter struct, then check T exists, is finite
% 4-by-1, and meets force/bound checks. Neither scale_dir=0 nor fallback=0 is
% a convergence certificate. fallback=1 is not a KKT optimality certificate.
% There is no previous-solution warm start. Method reference: Gouttefarde,
% Lamaury, Reichert and Bruckmann (2015), manuscript [29]; see METHOD_SOURCES.

    Tmin=WDM.param.lim_inf;
    Tmax=WDM.param.lim_sup;
    WDM.scale_dir=0;
    WDM.fallback=0;   % 1 = solution is the min-norm feasible vertex, not the KKT optimum
    WDM.param.M_eff=WDM.param.M;
    
    % Unit cable directions from the platform to the anchors.
    v = WDM.param.M_eff-repmat(P,1,size(WDM.param.M_eff,2));
    for i=1:size(WDM.param.M_eff,2)
        d=norm(v(:,i)); d=(d>WDM.param.TOLL)*d;
        if d==0 
            v(:,i)=[0 0]';
        else
            v(:,i)=v(:,i)/d;
        end
    end


    %% Determine of the Feasible Polygon Vertices

    % Compute N and t_p

    N = null(v);


    lookuptable = zeros(8,8,2);

    invV = get_pinv(v);
    t_p = invV * F;

    

    constraints = [Tmin - t_p; Tmax - t_p];
    vertices = [];

    i = 1;
    j = 2;
    b_list = [0 0 0];
    n_list = j;

    % find first vertex
    n_i = N(i,:);
    n_j = N(j,:);
    c_i = constraints(i);
    c_j = constraints(j);

    v_f = [n_i;n_j] \ [c_i;c_j];


    vertices = v_f;



    % find I
    I = find(constraints(1:4) - WDM.param.TOLL <= N * v_f & constraints(5:8) + WDM.param.TOLL >= N * v_f);

    iters = 0;
    while true
        % safety cap: on degenerate polygons the walk can cycle without ever
        % closing within TOLL; give up and return with no solution fields
        iters = iters + 1;
        if iters > 100
            return
        end
        % find n_i_orth
        n_i_orth = [N(i,2), -N(i,1)]';
        if N(j,:) * n_i_orth < 0
            if b_list(end-1) == 0
                n_i_orth = -n_i_orth;
            end
        else
            if b_list(end-1) == 1
                n_i_orth = -n_i_orth;
            end
        end

        % compute alpha_k for all k ~= i
        alpha_list = zeros(4,1);
        min_max_list = zeros(4,1);

        v_ij = [n_i;n_j] \ [c_i;c_j];

        for k = 1:4
            alpha = Inf;
            min_max = 0;
            n_v = N(k,:) * v_ij;
            n_k_n_i = N(k,:) * n_i_orth;

            if k == i
                % do nothing
            elseif abs(n_k_n_i) <= WDM.param.TOLL
                % do nothing
            elseif n_k_n_i > WDM.param.TOLL
                if n_v < constraints(k) - WDM.param.TOLL
                    alpha = (constraints(k) - n_v) / n_k_n_i;
                elseif constraints(k) - WDM.param.TOLL < n_v && n_v <= constraints(k+4) + WDM.param.TOLL
                    alpha = (constraints(k+4) - n_v) / n_k_n_i;
                    min_max = 1;
                else
                    % do nothing
                end
            elseif n_k_n_i < -WDM.param.TOLL
                if n_v <= constraints(k) - WDM.param.TOLL
                    % do nothing
                elseif constraints(k) - WDM.param.TOLL < n_v && n_v <= constraints(k+4) + WDM.param.TOLL
                    alpha = (constraints(k) - n_v) / n_k_n_i;
                else
                    alpha = (constraints(k+4) - n_v) / n_k_n_i;
                    min_max = 1;
                end
            else
                error("No Solution")
            end
            alpha_list(k) = alpha;
            min_max_list(k) = min_max;

        end

        [~,l] = min(alpha_list);
        b_list = [b_list, min_max_list(l)];
        n_list = [n_list, i];


        n_l = N(l,:);

        if b_list(end)
            c_l = constraints(l+4);
        else
            c_l = constraints(l);
        end

        v_li = [n_i;n_l] \ [c_i;c_l];


        lookuptable(i+4*b_list(end-1),l+4*b_list(end),:) = v_li;
        lookuptable(l+4*b_list(end),i+4*b_list(end-1),:) = v_li;

        if ~ismember(l,I)
            I = union(I,l);
            v_f = v_li;
            j = i;
            i = l;
            n_j = n_i;
            c_j = c_i;
            n_i = n_l;
            c_i = c_l;

        else
            
            if abs(v_li - v_f) < WDM.param.TOLL
                vertices = [vertices, v_li];
                break;
            else
                j = i;
                i = l;
                n_j = n_i;
                c_j = c_i;
                n_i = n_l;
                c_i = c_l;
            end
        end
        
        vertices = [vertices, v_li];
    end
    b_list = b_list(2:end-1);
    WDM.vertices = vertices;
    lambda = Inf;
    if length(I) == 4
        [n_list, b_list] = reduceListEfficient(n_list, b_list);

        for i = 1:length(n_list) - 1

            [a_i,~,~,~,index_i] = get_mu(i,n_list,b_list,N,t_p,Tmin,Tmax);
            [a_j,~,~,~,index_j] = get_mu(i+1,n_list,b_list,N,t_p,Tmin,Tmax);

            v_ij = [lookuptable(index_i, index_j,1); ...
                    lookuptable(index_i, index_j,2)];



            mu = get_inv([a_i', a_j']) * v_ij;

            if all(mu>0)
                lambda = v_ij;
                break
            end
        end

        if numel(lambda)==1 && lambda == Inf
            for i = 1:length(n_list)
                [a_i,~,~,mu_i,~] = get_mu(i,n_list,b_list,N,t_p,Tmin,Tmax);
                if mu_i >= 0
                    lambda = (mu_i * a_i)';
                    break
                end
            end
        end
    end

    if numel(lambda)==2
        % KKT optimum found (vertex or edge)
        WDM.T = N * lambda + t_p;
    else
        % ---- fallback: best (min-L2-norm) feasible traced vertex ----
        % The KKT vertex/edge search did not return a valid solution
        % (degenerate polygon: length(I)~=4, or no vertex/edge satisfied the
        % sign test). The feasible region was still traced, so we log the
        % lowest-norm vertex that respects the tension bounds instead of
        % discarding the point. This is an upper bound on the true L2-optimum.
        bestNorm = Inf; bestT = [];
        for kv = 1:size(vertices,2)
            Tk = N*vertices(:,kv) + t_p;
            if all(Tk >= Tmin - WDM.param.TOLL) && all(Tk <= Tmax + WDM.param.TOLL)
                nk = sqrt(Tk'*Tk);
                if nk < bestNorm
                    bestNorm = nk; bestT = Tk;
                end
            end
        end
        if ~isempty(bestT)
            WDM.T = bestT;
            WDM.fallback = 1;
        end
    end

    if isfield(WDM,'T') && numel(WDM.T)==4
        WDM.normT=sqrt(WDM.T'*WDM.T);
        WDM.errT=norm(v*WDM.T-F);
        if any(WDM.T < Tmin - WDM.param.TOLL | WDM.T > Tmax + WDM.param.TOLL)
            WDM.scale_dir = 1;
        end
    end



end

function [n_list, b_list] = reduceListEfficient(n_list, b_list)
    lastElement = n_list(end) + 4 * b_list(end);
    index = numel(n_list) - 1;
    found = false;
    
    while index > 0
        if n_list(index) + 4 * b_list(index) == lastElement
            found = true;
            break;
        end
        index = index - 1;
    end
    
    sortedIdx = index:length(n_list);
    
    % Extract the elements based on sorted indices
    n_list = n_list(sortedIdx);
    b_list = b_list(sortedIdx);
end

function [a,s,b,mu,index] = get_mu(i,n_list,b_list,N,t_p,Tmin,Tmax)
    index = n_list(i);
    if b_list(i) == 0
        b = Tmin(index);
        s = 1;
    else
        b = Tmax(index);
        s = -1;
    end
    a = s * N(index,:);
    mu = s * (b - t_p(index)) / (a * a');
end

function invV = get_pinv(v)
    V = v*v';
    a = V(1);
    b = V(2);
    c = V(3);
    d = V(4);

    detV = a*d - b*c;
    invV = v' * (1/detV) * [d, -b; -c, a];
end

function invV = get_inv(V)
    a = V(1,1);
    b = V(1,2);
    c = V(2,1);
    d = V(2,2);

    detV = a*d - b*c;
    invV = (1/detV) * [d, -b; -c, a];
end
