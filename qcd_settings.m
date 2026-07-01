data_prefix = "./data/qcd/";
figure_prefix = "./figure/qcd/qcd";

method_list = [...
    "FOM", ...
    "sFOM_s", ...
    "sGMRES_s", ...
    "adaFOM_t", ...
    "adaGMRES_t", ...
    ];

m = 150;
max_restarts = 300;

truncation_length_list = [2, 1, 0];

sketching_mat_type = "sparse sign";
sketching_size = 2 * m;
sketching_size_control = 2;
cond_tol = 1e6;

quad_tol = 1e-9;
stop_tol = 1e-7;
max_num_quad_points = 512;

number_thick = 5;

check_flag = 0;

fprintf("quad_tol: %.1e\n", quad_tol);
fprintf("stop_tol: %.1e\n", stop_tol);

fprintf("m: %d\n", m);
fprintf("max_restarts: %d\n", max_restarts);
fprintf("truncation_length: ");
for it = 1 : length(truncation_length_list)
    fprintf("%d ", truncation_length_list(it));
end
fprintf("\n");
fprintf("max_num_quad_points: %d\n", max_num_quad_points);

fprintf("sketching_mat_type: %s\n", sketching_mat_type);
fprintf("sketching_size: %d\n", sketching_size);
fprintf("sketcing_size_control: %d\n", sketching_size_control);
fprintf("cond_tol: %.1e\n", cond_tol);

fprintf("number_thick: %d\n", number_thick);

fprintf("data_prefix: %s\n", data_prefix);
fprintf("figure_prefix: %s\n", figure_prefix);

%% Load matrix.
load("./data/qcd_data/qcd_matrix_nonhermitian.mat");
N = size(Q, 1);
A = @(v) Q * (Q * v);
load("./data/qcd_data/qcd_nonhermitian_exact.mat");
c = zeros(N,1); c(1) = 1;
b = Q*c; normv = norm(b); b = b/norm(b);
f_ex = qcd_nonhermitian_exact / normv;

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
basic_param.exact = f_ex;
basic_param.stopping_accuracy = stop_tol;
basic_param.inner_product = @(a,b) b'*a;
basic_param.thick = [];
basic_param.min_decay = .95;
basic_param.waitbar = 0;
basic_param.reorth_number = 0;
basic_param.truncation_length = inf;
basic_param.verbose = 1;
basic_param.check = check_flag;

%% Print
caption_name = "Relative error, time, iteration number and matrix-vector product number" + ...
    " of the invsqrt function" + ...
    " for the QCD example.";
label_name = "qcd";
