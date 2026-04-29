clear all;
close all;

data_prefix = "./data/qcd_invsqrt/";
figure_prefix = "./figure/qcd_invsqrt/qcd_invsqrt";

method_list = [...
    "FOM", ...
    "sFOM_s", ...
    "sFOM_t"
    ];
truncation_length_list = [2, 1, 0];

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
fprintf("figure_prefix: %s\n", figure_prefix);

%% Load matrix.
load("./data/qcd/qcd_matrix_nonhermitian.mat");
load("./data/qcd/qcd_nonhermitian_exact.mat");
N = size(Q, 1);
A = @(v) Q * (Q * v);
c = zeros(N,1); c(1) = 1;
b = Q*c; normv = norm(b); b = b/norm(b);
f_ex = qcd_nonhermitian_exact / normv;

%% Plot results.
plot_figures(data_prefix, f_ex, ...
    method_list, truncation_length_list, ...
    figure_prefix);
