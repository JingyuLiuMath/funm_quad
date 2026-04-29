clear all;
close all;

data_prefix = "./data/network_exp/";
figure_prefix = "./figure/network_exp";

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
load('./data/network/wiki-Vote.mat');
load('./data/network/wiki-Vote-comp.mat');

A = -Problem.A; ee = -ee; N = size(A,1); 
b = ones(N,1);
f_ex = ex_expm;

%% Plot results.
plot_figures(data_prefix, f_ex, ...
    method_list, truncation_length_list, ...
    figure_prefix);
