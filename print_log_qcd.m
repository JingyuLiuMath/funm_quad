clear all;
close all;

truncation_length_list = [2, 1, 0];
method_list = [...
    "fix_FOM", ...
    "fix_sFOM_s", ...
    "ada_sFOM_t"
    ];

data_prefix = "./data/qcd_log/";

fprintf("truncation_length: ");
for it = 1 : length(truncation_length_list)
    fprintf("%d ", truncation_length_list(it));
end
fprintf("\n");
fprintf("method: ")
for it = 1 : length(method_list)
    fprintf("%s ", method_list(it));
end
fprintf("\n");
fprintf("data_prefix: %s\n", data_prefix);


%% Load matrix.
load("./data/qcd/qcd_matrix_nonhermitian.mat");
N = size(Q, 1);
A = @(v) Q * (Q * v);
c = zeros(N,1); c(1) = 1;
b = Q*c; normv = norm(b); b = b/norm(b);

addpath('funm_quad')
exact_param.function = 'log';
exact_param.restart_length = 300; 
exact_param.max_restarts = 15;
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

%% print methods.
caption_name = "Relative error, time, iteration number and oracle number" + ...
    " of the log function" + ...
    " for the QCD example";
label_name = "log_qcd";

print_table(data_prefix, f_ex, ...
    method_list, ...
    truncation_length_list, ...
    caption_name, label_name);

