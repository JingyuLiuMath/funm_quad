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

file_prefix = "./data/frac_laplacian_invsqrt/";

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

%% print methods.
caption_name = "Relative error, time and oracle number" + ...
    " of the invsqrt function" + ...
    " for the fractional Laplacian example";
label_name = "tab:invsqrt_frac_laplacian";
fprintf("\n");
fprintf("\\begin{table}[tbhp]\n")
fprintf("\\centering\n")
fprintf("\\begin{tabular}{cccc}\n")
fprintf("\\toprule\n")
fprintf("Method & error & time & oracle number \\\\ \n")
num_truncate_len = length(truncation_length_list);
num_method = length(method_list);
for it_method = 1 : num_method
    method_name = method_list(it_method);
    if endsWith(method_name, "_t")
        for it_trunc = 1 : num_truncate_len
            truncation_length = truncation_length_list(it_trunc);
            save_name = method_name + "_" + string(truncation_length);
            load(file_prefix + save_name  + ".mat");
            print_result(save_name, result, f_ex);
        end
    else
        truncation_length = inf;
        save_name = method_name;
        load(file_prefix + save_name  + ".mat");
        print_result(save_name, result, f_ex);
        
    end
end
fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\caption{%s}\n", caption_name);
fprintf("\\label{tab:%s}\n", label_name);
fprintf("\\end{table}\n");
fprintf("\n");

function print_result(save_name, result, f_ex)

rel_err = norm(f_ex - result.f) / norm(f_ex);
print_name = replace(save_name, "_", "-");
fprintf("\\midrule\n")
fprintf("%s ", print_name );
fprintf("& %.1e ", rel_err);
fprintf("& %.1e ", result.time);
fprintf("& %d ", result.num_oracle);
fprintf("\\\\ \n");

end