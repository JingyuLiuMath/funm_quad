clear;
close all;
rng(2026);
maxNumCompThreads(1);
warning off;

quad_tol = 1e-7;
stop_tol = 1e-8;

m = 300;
max_restarts = 15;
ada_max_restarts = 200;
truncation_length = [2, 1, 0];
max_num_quad_points = 8192;

sketching_mat_type = "sparse sign";
sketching_size = 2 * m;
ada_sketching_size_control = 2;
cond_tol = 1e4;

method_list = ["benchmark", ...
    "fix FOM-s", ...
    "ada FOM-t last-orth", "ada FOM-t last-sorth", "ada FOM-t whitening", ...
    "ada FOM-s last-orth", "ada FOM-s last-sorth", "ada FOM-s whitening", ...
    ];

save_flag = 1;
check_flag = 0;

fprintf("quad_tol: %.1e\n", quad_tol);
fprintf("stop_tol: %.1e\n", stop_tol);

fprintf("m: %d\n", m);
fprintf("max_restarts: %d\n", max_restarts);
fprintf("ada_max_restarts: %d\n", ada_max_restarts);
fprintf("truncation_length: %s\n", format_numeric_vector(truncation_length));
fprintf("max_num_quad_points: %d\n", max_num_quad_points);

fprintf("sketching_mat_type: %s\n", sketching_mat_type);
fprintf("sketching_size: %d\n", sketching_size);
fprintf("ada_sketcing_size_control: %d\n", ada_sketching_size_control);
fprintf("cond_tol: %.1e\n", cond_tol);

%% Load matrix.
load('./data/frac_laplacian/gnutella_comp.mat');
A = L;
N = size(A, 1);
b = A * b;
f_ex = ex_gnutella;

%% choose parameters for the FUNM_QUAD restart algorithm
% jingyu: tol and stopping_accruacy are modified
addpath('funm_quad')
basic_param.function = 'invSqrt';
basic_param.restart_length = m;          % each restart cycle consists of 70 Arnoldi iterations
basic_param.max_restarts = max_restarts;            % perform at most 15 restart cycles
basic_param.tol = quad_tol;                   % tolerance for quadrature rule
basic_param.transformation_parameter = 1;     % parameter for the integral transformation
basic_param.hermitian = 0;                % the matrix A is Hermitian
basic_param.V_full = 0;                   % set 1 if you need Krylov basis
basic_param.H_full = 0;                   % do not store all Hessenberg matrices
basic_param.exact = [];
basic_param.stopping_accuracy = stop_tol;     % stopping accuracy
basic_param.inner_product = @(a,b) b'*a;  % use standard Euclidean inner product
basic_param.thick = [];                   % no implicit deflation is performed
basic_param.min_decay = .95;              % we desire linear error reduction of rate < .95 
basic_param.waitbar = 0;                  % show waitbar 
basic_param.reorth_number = 0;              % #reorthogonalizations
basic_param.truncation_length = inf;      % truncation length for Arnoldi 
basic_param.verbose = 1;                  % print information about progress of algorithm
basic_param.check = check_flag;

add_param = construct_ada_param(...
    truncation_length, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    ada_max_restarts, ada_sketching_size_control, cond_tol);

%% test methods
result_list = run_methods(A, b, f_ex, method_list, basic_param, add_param);

%% save data
file_name = "./data/frac_laplacian/frac_laplacian_.mat";
if save_flag
    save(file_name, "result_list");
end

%% print table
fprintf("\n\n");
print_table(result_list);

%% plot
exmaple_name = "frac_laplacian";
save_path = "./figure/" + exmaple_name + "/";
plot_figures(result_list, exmaple_name, save_path, save_flag);

function text = format_numeric_vector(values)
text = "[" + strjoin(string(values), ", ") + "]";
end
