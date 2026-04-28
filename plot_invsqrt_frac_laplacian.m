clear all;
close all;

method_list = [...
    "FOM", ...
    "sFOM_s", ...
    "sFOM_t"
    ];
truncation_length_list = [2, 1, 0];

data_prefix = "./data/frac_laplacian_invsqrt/";
figure_prefix = "./figure/frac_laplacian_invsqrt/frac_laplacian_invsqrt";

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
load('./data/frac_laplacian/gnutella_comp.mat');
f_ex = ex_gnutella;

%% Plot results.
plot_figures(data_prefix, f_ex, ...
    method_list, truncation_length_list, ...
    figure_prefix);
