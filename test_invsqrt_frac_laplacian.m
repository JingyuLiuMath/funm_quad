clear all;
close all;
rng(2026);
maxNumCompThreads(1);
warning off;

method_list = [...
    "FOM", ...
    "sFOM_s", ...
    "sFOM_t"
    ];

m = 300;
max_restarts = 300;

truncation_length_list = [2, 1, 0];

sketching_mat_type = "sparse sign";
sketching_size = 2 * m;
sketching_size_control = 2;
cond_tol = 1e4;

quad_tol = 1e-7;
stop_tol = 1e-8;
max_num_quad_points = 8192;

check_flag = 0;

data_prefix = "./data/frac_laplacian_invsqrt/";

fprintf("quad_tol: %.1e\n", quad_tol);
fprintf("stop_tol: %.1e\n", stop_tol);

fprintf("m: %d\n", m);
fprintf("max_restarts: %d\n", max_restarts);
fprintf("ada_max_restarts: %d\n", max_restarts);
fprintf("truncation_length: ");
for it = 1 : length(truncation_length_list)
    fprintf("%d ", truncation_length_list(it));
end
fprintf("\n");
fprintf("max_num_quad_points: %d\n", max_num_quad_points);

fprintf("sketching_mat_type: %s\n", sketching_mat_type);
fprintf("sketching_size: %d\n", sketching_size);
fprintf("ada_sketcing_size_control: %d\n", sketching_size_control);
fprintf("cond_tol: %.1e\n", cond_tol);

fprintf("file_prefix: %s\n", data_prefix);

%% Load matrix.
load('./data/frac_laplacian/gnutella_comp.mat');
A = L;
N = size(A, 1);
b = A * b;
f_ex = ex_gnutella;

%% choose parameters for the FUNM_QUAD restart algorithm
addpath('funm_quad')
basic_param.function = 'invSqrt';
basic_param.restart_length = m;
basic_param.max_restarts = max_restarts;
basic_param.tol = quad_tol;
basic_param.transformation_parameter = 1;
basic_param.hermitian = 0;
basic_param.V_full = 0;
basic_param.H_full = 0;
basic_param.exact = [];
basic_param.stopping_accuracy = stop_tol;
basic_param.inner_product = @(a,b) b'*a;
basic_param.thick = [];
basic_param.min_decay = .95;
basic_param.waitbar = 0;
basic_param.reorth_number = 0;
basic_param.truncation_length = inf;
basic_param.verbose = 1;
basic_param.check = check_flag;

%% Test methods
run_methods(data_prefix, ...
    method_list, ...
    truncation_length_list, ...
    basic_param, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    sketching_size_control, cond_tol, ...
    A, b);
