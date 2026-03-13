clear;
close all;
rng(2026);
maxNumCompThreads(1);
warning off;

quad_tol = 1e-7;
stop_tol = 1e-8;

m = 100;
max_restarts = 15;
ada_max_restarts = 200;
truncation_length = 2;
max_num_quad_points = 8192;

sketching_mat_type = "sparse sign";
sketching_size = 2 * m;
ada_sketching_size_control = 2;
cond_tol = 1e4;

method_list = ["benchmark", ...
    "fix FOM-t whitening", ...
    "fix FOM-s", ...
    "ada FOM-t last-orth", "ada FOM-t last-sorth", "ada FOM-t whitening", ...
    "ada FOM-s whitening", "ada FOM-st whitening"];

fprintf("quad_tol: %.4e\n", quad_tol);
fprintf("stop_tol: %.4e\n", stop_tol);

fprintf("m: %d\n", m);
fprintf("max_restarts: %d\n", max_restarts);
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
basic_param.max_restarts = max_restarts;            % perform at most 15 restart cycles
basic_param.tol = quad_tol;                   % tolerance for quadrature rule
basic_param.hermitian = 0;                % the matrix A is Hermitian
basic_param.V_full = 0;                   % set 1 if you need Krylov basis
basic_param.H_full = 0;                   % do not store all Hessenberg matrices
basic_param.exact = [];
basic_param.stopping_accuracy = stop_tol;     % stopping accuracy
basic_param.inner_product = @(a,b) b'*a;  % use standard Euclidean inner product
basic_param.thick = [];                   % no implicit deflation is performed
basic_param.min_decay = .95;              % we desire linear error reduction of rate < .95
basic_param.waitbar = 0;                  % show waitbar
basic_param.reorth_number = 0;            % #reorthogonalizations
basic_param.truncation_length = inf;      % truncation length for Arnoldi
basic_param.verbose = 1;                  % print information about progress of algorithm

add_param = construct_ada_param(...
    truncation_length, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    ada_max_restarts, ada_sketching_size_control, cond_tol);

%% test methods
result_list = run_methods(A, b, f_ex, method_list, basic_param, add_param);

%% save data
file_name = "./data/wikivote/wikivote_trunc_" + string(truncation_length) + ".mat";
save(file_name, "result_list");

%% print table
fprintf("\n\n");

fprintf("t = %d\n", truncation_length);

print_table(result_list);

%% plot 
save_path = "./figure/wikivote/";
plot_figures(result_list, save_path, truncation_length);
