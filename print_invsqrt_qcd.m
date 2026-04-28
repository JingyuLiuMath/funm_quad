clear all;
close all;

truncation_length_list = [2, 1, 0];
method_list = [...
    "fix_FOM", ...
    "fix_sFOM_s", ...
    "ada_sFOM_t"
    ];

data_prefix = "./data/qcd_invsqrt/";

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
load("./data/qcd/qcd_nonhermitian_exact.mat");
c = zeros(size(Q,1),1); c(1) = 1;
b = Q*c; normv = norm(b);
f_ex = qcd_nonhermitian_exact / normv;

%% print methods.
caption_name = "Relative error, time, iteration number and oracle number" + ...
    " of the invsqrt function" + ...
    " for the QCD example";
label_name = "invsqrt_qcd";

print_table(data_prefix, f_ex, ...
    method_list, ...
    truncation_length_list, ...
    caption_name, label_name);

