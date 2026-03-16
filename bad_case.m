clear;
close all;
rng(2026);
maxNumCompThreads(1);

% m = 50;
% m_max = 100;
m = 100;
m_max = 200;
ada_max_restarts = 15;
truncation_length = 0;
max_num_quad_points = 8192;

sketching_mat_type = "Clarkson-Woodruff";
sketching_size = 2 * m;
ada_sketching_size_control = 2;
cond_tol = 1e4;

method_name = "ada FOM-t last-sorth";

check_flag = 1;

fprintf("m: %d\n", m);
fprintf("m_max: %d\n", m_max);
fprintf("ada_max_restarts: %d\n", ada_max_restarts);
fprintf("truncation_length: %d\n", truncation_length);
fprintf("max_num_quad_points: %d\n", max_num_quad_points);

fprintf("sketching_mat_type: %s\n", sketching_mat_type);
fprintf("sketching_size: %d\n", sketching_size);
fprintf("ada_sketcing_size_control: %d\n", ada_sketching_size_control);
fprintf("cond_tol: %d\n", cond_tol);

%% Load matrix.
load('data/wikivote/wiki-Vote.mat');
load('data/wikivote/wiki-Vote-comp.mat');
A = -Problem.A; ee = -ee; N = size(A,1); 
b = ones(N,1);
f_ex = ex_expm;

%% choose parameters for the FUNM_QUAD restart algorithm
% jingyu: tol and stopping_accruacy are modified
addpath('funm_quad')
basic_param.function = 'exp';
basic_param.restart_length = m;          % each restart cycle consists of 70 Arnoldi iterations
basic_param.max_restarts = 15;            % perform at most 15 restart cycles
basic_param.tol = 1e-10;                   % tolerance for quadrature rule
basic_param.hermitian = 0;                % the matrix A is Hermitian
basic_param.V_full = 0;                   % set 1 if you need Krylov basis
basic_param.H_full = 0;                   % do not store all Hessenberg matrices
basic_param.exact = [];
basic_param.stopping_accuracy = 1e-14;     % stopping accuracy
basic_param.inner_product = @(a,b) b'*a;  % use standard Euclidean inner product
basic_param.thick = [];                   % no implicit deflation is performed
basic_param.min_decay = .95;              % we desire linear error reduction of rate < .95
basic_param.waitbar = 0;                  % show waitbar
basic_param.reorth_number = 0;            % #reorthogonalizations
basic_param.truncation_length = inf;      % truncation length for Arnoldi
basic_param.verbose = 2;                  % print information about progress of algorithm
basic_param.check = check_flag;

add_param = construct_ada_param(...
    truncation_length, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    ada_max_restarts, ada_sketching_size_control, cond_tol);

%% afom-t
rng(2029)

result = run_single_method(A, b, f_ex, method_name, basic_param, add_param);

f = result{1}.f;
out = result{1}.out;
t = result{1}.time;

num_it = size(out.appr, 2);
rel_err = norm(f_ex - f) / norm(f_ex);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t);
at_rel_err0 = norm(out.appr(:, 1) - f_ex) / norm(f_ex);
fprintf("initial err: %e\n", at_rel_err0);
fprintf("\n\n");
iter_err = vecnorm(f_ex - out.appr) / norm(f_ex);
fprintf("rel err in each iter\n");
for j = 1 : length(iter_err)
    fprintf("%.4e ", iter_err(j));
end
fprintf("\n");

fprintf("num quad points in each iter\n");
for j = 1 : length(iter_err)
    fprintf("%d ", out.num_quadpoints(j));
end
fprintf("\n");


fprintf("subspace dim in each iter\n");
for j = 1 : length(iter_err)
    fprintf("%d ", out.dim(j));
end
fprintf("\n");
