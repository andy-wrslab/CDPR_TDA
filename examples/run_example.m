function result = run_example(mode)
%RUN_EXAMPLE Display and check one planar four-cable tension distribution.
%   RESULT = RUN_EXAMPLE runs the proposed Direct Method (DM). It requires
%   only MATLAB. RESULT = RUN_EXAMPLE('all') also runs the five comparison
%   methods and WEC; bounded QP requires Optimization Toolbox.
%
%   MODE is 'direct' (default) or 'all'. The example creates a known feasible
%   force at a strictly interior position. Anchor coordinates and position
%   are in mm; requested/achieved forces and four cable tensions are in N.
%   Anchor columns are counterclockwise; tension rows use that same order.
%
%   RESULT contains position_mm (2x1), anchors_mm (2x4), requested_force_N
%   and wec_target_N (2x1), lower/upper_bounds_N (4x1), and a methods struct
%   array. Each method has tensions_N (4x1), achieved_force_N (2x1), scalar
%   equilibrium_error_N and target_error_N (Euclidean norms), bounds_ok,
%   finite_output, native_exitflag (NaN when absent), native_ok and passed.
%   native_ok means a positive exitflag when one exists; it is not an
%   independent proof of force reproduction. passed combines finite output,
%   bounds (1e-9 N allowance), native status, and raw-force error <=1e-8 N.
%   native_output preserves the original returned struct for inspection.
%   This example does not measure execution times or modify solver settings.

if nargin < 1, mode = 'direct'; end
mode = validatestring(mode, {'direct','all'}, mfilename, 'mode');
if strcmp(mode,'all') && (~license('test','Optimization_Toolbox') || ...
        exist('quadprog','file') ~= 2)
    error('run_example:MissingToolbox', ...
        'The all-method example requires Optimization Toolbox. Use run_example for DM only.');
end
root = fileparts(fileparts(mfilename('fullpath')));
original_path = path;
path_cleanup = onCleanup(@() path(original_path)); %#ok<NASGU>
addpath(fullfile(root,'src'));

% Fixed CCW cable order, starting at the bottom-left anchor.
parameters.param.M = [-315 315 315 -315; -315 -315 315 315];
parameters.param.lim_inf = 3.8*ones(4,1);
parameters.param.lim_sup = 25*ones(4,1);
parameters.param.TOLL = 1e-14;
parameters.param.TOLL_WEC = 1e-9;
parameters.verbose = false;
parameters.acta = struct('tol',1e-10);
parameters.xacta = struct('tol',1e-9,'force_tol',1e-8);
position = [80;45];
directions = parameters.param.M-position;
directions = directions./sqrt(sum(directions.^2,1));
% Construct a feasible request; the method need not return these seed tensions.
requested_force = directions*[12;13;14;15];
target_force = requested_force;
wec_flag = NaN;

names = {'DM'};
functions = {@wt4_2024_minmax_v4};
if strcmp(mode,'all')
    [target_force,wec_flag] = WEC_v5(parameters,position,requested_force);
    assert(isequal(size(target_force),[2 1]) && all(isfinite(target_force)), ...
        'run_example:UnavailableTarget','WEC did not return a finite force for the feasible example.');
    names = {'DM','ICFM','Bounded QP','VTDA-L2','ACTA','X-ACTA'};
    functions = {@wt4_2024_minmax_v4,@wt_pott_v2,@wt_qp_bounded_v11, ...
        @wt_gouttefarde_v2,@wt_acta,@wt_xacta};
end
result = struct('position_mm',position,'anchors_mm',parameters.param.M, ...
    'requested_force_N',requested_force,'wec_target_N',target_force, ...
    'lower_bounds_N',parameters.param.lim_inf, ...
    'upper_bounds_N',parameters.param.lim_sup, ...
    'wec_native_flag',wec_flag,'methods',struct([]));

fprintf('\nTension distribution example\n');
fprintf('Position [x y]:         [%g %g] mm\n',position);
fprintf('Requested force [x y]:  [% .9f % .9f] N\n',requested_force);
fprintf('Tension bounds:        3.8 <= T_i <= 25 N\n');
if strcmp(mode,'all')
    fprintf('Shared WEC target:     [% .9f % .9f] N\n',target_force);
end
for k = 1:numel(functions)
    input_force = target_force;
    if strcmp(names{k},'X-ACTA'), input_force = requested_force; end
    % Fresh parameters prevent a prior returned tension vector being reused.
    output = functions{k}(parameters,position,input_force);
    finite_output = isfield(output,'T') && isequal(size(output.T),[4 1]) ...
        && all(isfinite(output.T));
    tensions = NaN(4,1); achieved = NaN(2,1);
    bound_error = Inf; equilibrium_error = Inf; target_error = Inf;
    if finite_output
        tensions = output.T;
        achieved = directions*tensions;
        bound_error = max([0; parameters.param.lim_inf-tensions; ...
            tensions-parameters.param.lim_sup]);
        equilibrium_error = norm(achieved-requested_force,2);
        target_error = norm(achieved-input_force,2);
    end
    native_flag = NaN; native_ok = true;
    if isfield(output,'exitflag')
        native_flag = output.exitflag;
        native_ok = isscalar(native_flag) && isfinite(native_flag) && native_flag > 0;
    end
    bounds_ok = finite_output && bound_error <= 1e-9;
    passed = bounds_ok && native_ok && equilibrium_error <= 1e-8;
    item = struct('name',names{k},'tensions_N',tensions, ...
        'achieved_force_N',achieved,'equilibrium_error_N',equilibrium_error, ...
        'target_error_N',target_error,'bound_violation_N',bound_error, ...
        'bounds_ok',bounds_ok,'finite_output',finite_output, ...
        'native_exitflag',native_flag,'native_ok',native_ok, ...
        'passed',passed,'native_output',output);
    if k == 1
        result.methods = item;
    else
        result.methods(k) = item;
    end
    fprintf('\n%s\n',names{k});
    fprintf('  Tensions [T1 T2 T3 T4]: [% .9f % .9f % .9f % .9f] N\n',tensions);
    fprintf('  Achieved force [x y]:   [% .9f % .9f] N\n',achieved);
    fprintf('  Equilibrium error:     %.3e N\n',equilibrium_error);
    fprintf('  Tension bounds satisfied: %s\n',yes_no(bounds_ok));
    if isfield(output,'exitflag'), fprintf('  Native exit flag:      %g\n',native_flag); end
    fprintf('  All example checks passed: %s\n',yes_no(passed));
end
fprintf('\nThis example produces no benchmark timing observations.\n');
end

function value = yes_no(condition)
if condition, value = 'yes'; else, value = 'no'; end
end
