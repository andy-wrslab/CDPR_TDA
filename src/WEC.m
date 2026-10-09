function [F_new, is_scaled] = WEC(WDM,P,F)
% WEC Adjust a requested force using wrench exertion capability.
%   [F_NEW,IS_SCALED] = WEC(WDM,P,F) maps the 16 tension-box vertices
%   into force space and constructs their convex hull. For nonzero F it
%   intersects the SIGNED force line with that polygon and adjusts the
%   force by the inward TOLL_WEC margin when needed.
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
% Also requires WDM.param.TOLL_WEC, an inward force margin [N] (1e-9 in the
% reviewed setup). Equal cable bounds are required: only their first entries
% are used to build the four-cable box. No tensions are returned by WEC.
%
% OUTPUTS
%   F_NEW       2-by-1 adjusted force [N], or NaN(2,1) when the nonzero
%               requested force line does not intersect the attainable set.
%   IS_SCALED   Legacy logical flag: false only for the unchanged nonzero
%               branch; true on other branches, including unchanged zero.
% Compute norm(F_NEW-F,2) separately to determine actual force change. Signed
% line adjustment can reduce, increase, or reverse F; it is not constrained
% to a positive scaling factor. An infeasible zero command receives a nonzero
% fallback. Degenerate hull/geometry errors may propagate. A finite result
% still needs an independent attainable-force check.
%
% With verbose=true, a diagnostic figure is drawn and execution pauses; use
% false for numerical calls. Geometry/hull are rebuilt, with no lookup cache.
% Method context: Boschetti, Passarini, Trevisani and Zanotto (2018), A fast
% algorithm for wrench exertion capability computation, manuscript [34].
% See METHOD_SOURCES.txt; the captured file does not establish author-code origin.

is_scaled = true;
v = WDM.param.M-repmat(P,1,size(WDM.param.M,2));
Tmin=WDM.param.lim_inf(1);
Tmax=WDM.param.lim_sup(1);

for i=1:size(WDM.param.M,2)
    d=norm(v(:,i)); d=(d>WDM.param.TOLL)*d;
    if d==0
        v(:,i)=[0 0]';
    else
        v(:,i)=v(:,i)/d;
    end
end

nr=norm(F,2); nr=(nr>WDM.param.TOLL)*nr;
if nr ==0
    F_new=[0;0];
    isnull_F=true;
    vr=v;
    Rrb=eye(2); 
    %WDM.scale_dir=0;
else
    isnull_F=false;
    ur=F/nr;
    Rrb=[ur(1) -ur(2); ur(2) ur(1)];
    vr=Rrb'*v; % rotated struct matrix
end

%% define orthotope

% V_char = dec2bin([0:(2^4-1)],4)';
% V=nan(size(V_char));
% for j=1:size(V_char,2)
%     for i=1:size(V_char,1)
%         eval(['V(i,j)=[' V_char(i,j) ']'';']);
%     end
% end
% idx=logical(V);
% V(idx)=Tmax;
% V(~idx)=Tmin;

% this should save some time (define vertices explicitly)
V =[
    Tmin     Tmin     Tmin     Tmin     Tmin     Tmin     Tmin     Tmin     Tmax     Tmax     Tmax     Tmax     Tmax     Tmax     Tmax     Tmax
    Tmin     Tmin     Tmin     Tmin     Tmax     Tmax     Tmax     Tmax     Tmin     Tmin     Tmin     Tmin     Tmax     Tmax     Tmax     Tmax
    Tmin     Tmin     Tmax     Tmax     Tmin     Tmin     Tmax     Tmax     Tmin     Tmin     Tmax     Tmax     Tmin     Tmin     Tmax     Tmax
    Tmin     Tmax     Tmin     Tmax     Tmin     Tmax     Tmin     Tmax     Tmin     Tmax     Tmin     Tmax     Tmin     Tmax     Tmin     Tmax
    ];


% project orthotope to get candidate vertices of polytope
U_all=vr*V;

% extract actual vertices
[idx,~] = convhull(U_all'); %indices are sorted CCW
U=U_all(:,idx(1:end-1));

% get normals and characteristic vector delta
n=nan(2,size(U,2)-1);
delta=n(1,:)';
t=n;
for j=1:size(U,2)
    P1=U(:,j);
    P2=U(:,mod(j,size(U,2))+1);
    t(:,j)=P2-P1; t(:,j)=t(:,j)/norm(t(:,j),2);
    n(:,j)=[0 1; -1 0]*t(:,j);
    delta(j)=n(:,j)'*P1;
end
n=n'; % transpose
t=t';

if ~(isnull_F)
    %% nontrivial F

    c1=sum(U(2,:)>0)==size(U,2);
    c2=sum(U(2,:)<0)==size(U,2);
    % does the desired direction intersect the polygon?
    if c1 || c2
        warning('Cannot exert any F in the desired direction.');
        F_new=nan(size(F));
        return
    end
    %
    % compute subsets of edges:
    %n(:,1)=(abs(n(:,1))>WDM.param.TOLL).*n(:,1);
    idx_P=n(:,1)>0;
    idx_Q=n(:,1)<0;
    %idx_S=n(:,1)==0;

    lambda_wmax=min(delta(idx_P)./n(idx_P,1));
    lambda_wmin=max(delta(idx_Q)./n(idx_Q,1));

    % compute max and min Force projections along the desired line of action, add tolerances
    Fmax_proj=max([lambda_wmin,lambda_wmax]);
    Fmin_proj=min([lambda_wmin,lambda_wmax]);
    Fmax_safe=Fmax_proj-WDM.param.TOLL_WEC;
    Fmin_safe=Fmin_proj+WDM.param.TOLL_WEC;
    % correct force F, if necessary
    if (nr>=Fmin_safe) && (nr<=Fmax_safe)
        % F is feasible
        F_new=F;
        is_scaled = false;
        %disp('not changed')
        %disp('case A1')
    elseif (nr<Fmin_safe)
        F_new=Fmin_safe*Rrb*[1;0];
        % WDM.verbose=true;
        %disp('changed to Fmin')
        %disp('case A2')
    else
        F_new=Fmax_safe*Rrb*[1;0];
        % WDM.verbose=true;
        %disp('changed to Fmax')
        %disp('case A3')
    end
elseif (isnull_F) && (sum(delta>=0)==length(delta))
    %% origin inside, null force - do nothing
    %disp('case B')
    % disp(char(13));
else
    %% origin outside, null force – pick closest point
    [dist, closestPoint] = closestPointToOrigin(U');
    F_new=(dist+WDM.param.TOLL_WEC)*(closestPoint'/dist); %safe point
    %disp('case C')
    % disp(F)
    % disp(F_new)
    % disp(char(13));
    % WDM.verbose=true;
end


%% debug plot
if WDM.verbose
    myfig=figure;
    h1=plot(U_all(1,:),U_all(2,:),'ok');hold on;
    h2=patch('XData',U(1,:),'YData',U(2,:),'FaceColor', [211, 211, 211]/255,'EdgeColor','none','FaceAlpha',0.5);
    for j=1:size(U,2)
        ht(j)=text(U(1,j),U(2,j),num2str(j),'FontName','Arial','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','top');
        P1=U(:,j); P2=U(:,mod(j,size(U,2))+1); Pm=.5*(P1+P2);
        hn(j)=plot([Pm(1),Pm(1)+n(j,1)],[Pm(2),Pm(2)+n(j,2)],'-g');
        ht(j)=plot([Pm(1),Pm(1)+t(j,1)],[Pm(2),Pm(2)+t(j,2)],'-r');
    end
    h5=plot([0 nr],[0, 0],'-g');
    if ~isnull_F
        h3=plot([0 Fmax_proj],[0, 0],'-r');
        h4=plot([0 Fmin_proj],[0, 0],'-b');
    end
    h6=plot([1;0]'*Rrb'*F_new,[0;1]'*Rrb'*F_new,'*g');
    box on; grid on; axis equal;
    pause;
    close(myfig);
end
end


function [minDistance, closestPoint] = closestPointToOrigin(vertices)
    % vertices is an Nx2 matrix where N is the number of vertices
    % Each row represents the [x, y] coordinates of a vertex

    numVertices = size(vertices, 1);
    minDistance = inf; % Initialize minimum distance with infinity
    closestPoint = [0, 0]; % Initialize closest point

    for i = 1:numVertices
        % Calculate distance from the origin to the current vertex
        vertexDist = norm(vertices(i, :));
        if vertexDist < minDistance
            minDistance = vertexDist;
            closestPoint = vertices(i, :);
        end
        
        % Calculate distance from the origin to the edge segment
        if i < numVertices
            nextVertex = vertices(i + 1, :);
        else
            nextVertex = vertices(1, :); % Loop back to the first vertex for the last edge
        end
        
        edgeVec = nextVertex - vertices(i, :);
        pointVec = -vertices(i, :); % Vector from origin to the current vertex
        
        % Project point vector onto edge vector
        projLength = dot(pointVec, edgeVec) / norm(edgeVec)^2;
        if projLength > 0 && projLength < 1
            % The perpendicular projection is within the segment
            projection = vertices(i, :) + projLength * edgeVec;
            distToProjection = norm(projection);
            if distToProjection < minDistance
                minDistance = distToProjection;
                closestPoint = projection;
            end
        end
    end
end
