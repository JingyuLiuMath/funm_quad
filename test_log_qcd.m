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
truncation_length_list = [2, 1, 0];
max_num_quad_points = 8192;

sketching_mat_type = "sparse sign";
sketching_size = 2 * m;
ada_sketching_size_control = 2;
cond_tol = 1e4;

method_list = [...
    "fix_FOM", ...
    "fix_sFOM_s", ...
    "ada_sFOM_t"
    ];

check_flag = 0;

file_prefix = "./data/qcd_log/";

fprintf("quad_tol: %.1e\n", quad_tol);
fprintf("stop_tol: %.1e\n", stop_tol);

fprintf("m: %d\n", m);
fprintf("max_restarts: %d\n", max_restarts);
fprintf("ada_max_restarts: %d\n", ada_max_restarts);
fprintf("truncation_length: ");
for it = 1 : length(truncation_length_list)
    fprintf("%d ", truncation_length_list(it));
end
fprintf("\n");
fprintf("max_num_quad_points: %d\n", max_num_quad_points);

fprintf("sketching_mat_type: %s\n", sketching_mat_type);
fprintf("sketching_size: %d\n", sketching_size);
fprintf("ada_sketcing_size_control: %d\n", ada_sketching_size_control);
fprintf("cond_tol: %.1e\n", cond_tol);

fprintf("file_prefix: %s\n", file_prefix);

%% Load matrix.
load("./data/qcd/qcd_matrix_nonhermitian.mat");
N = size(Q, 1);
A = @(v) Q * (Q * v);
load("./data/qcd/qcd_nonhermitian_exact.mat");
c = zeros(N,1); c(1) = 1;
b = Q*c; normv = norm(b); b = b/norm(b);

exact_param.function = 'log';
exact_param.restart_length = m; 
exact_param.max_restarts = max_restarts;
exact_param.tol = 1e-15;
exact_param.transformation_parameter = 1;
exact_param.hermitian = 0;
exact_param.V_full = 0;
exact_param.H_full = 0;
exact_param.exact = [];
exact_param.stopping_accuracy = 1e-15;
exact_param.inner_product = @(a,b) b'*a;
exact_param.thick = [];
exact_param.min_decay = .95;
exact_param.waitbar = 0;
exact_param.reorth_number = 0;
exact_param.truncation_length = inf;
exact_param.verbose = 1;
[f_ex, ~] = funm_quad(A, b, exact_param);


%% choose parameters for the FUNM_QUAD restart algorithm
addpath('funm_quad')
basic_param.function = 'log';
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

%% test methods
num_truncate_len = length(truncation_length_list);
num_method = length(method_list);
for it_method = 1 : num_method
    method_name = method_list(it_method);
    if endsWith(method_name, "_t")
        for it_trunc = 1 : num_truncate_len
            truncation_length = truncation_length_list(it_trunc);
            
            save_name = method_name + "_" + string(truncation_length);
            fprintf("\n\n");
            fprintf("%s\n", save_name);

            add_param = construct_ada_param(...
                truncation_length, ...
                max_num_quad_points, ...
                sketching_mat_type, sketching_size, ...
                ada_max_restarts, ada_sketching_size_control, cond_tol);

            result = run_single_method(A, b, method_name, basic_param, add_param);
            result.save_name = save_name;

            save(file_prefix + save_name  + ".mat", "result");
        end
    else
        truncation_length = inf;

        save_name = method_name;
        fprintf("\n\n");
        fprintf("%s\n", save_name);
        
        add_param = construct_ada_param(...
            truncation_length, ...
            max_num_quad_points, ...
            sketching_mat_type, sketching_size, ...
            ada_max_restarts, ada_sketching_size_control, cond_tol);

        result = run_single_method(A, b, method_name, basic_param, add_param);
        result.save_name = save_name;
        
        save(file_prefix + save_name  + ".mat", "result");
    end
end

function text = format_numeric_vector(values)
text = "[" + strjoin(string(values), ", ") + "]";
end
