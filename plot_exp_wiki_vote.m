clear all;
close all;

truncation_length_list = [2, 1, 0];
method_list = [...
    "fix_FOM", ...
    "fix_sFOM_s", ...
    "ada_sFOM_t"
    ];

data_prefix = "./data/wiki_vote_exp/";
figure_prefix = "./figure/wiki_vote_exp/wiki_vote_exp";

results_list = [];

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
load('./data/wiki_vote/wiki-Vote.mat');
load('./data/wiki_vote/wiki-Vote-comp.mat');

A = -Problem.A; ee = -ee; N = size(A,1); 
b = ones(N,1);
f_ex = ex_expm;

%% Plot results.
plot_figures(data_prefix, f_ex, ...
    method_list, truncation_length_list, ...
    figure_prefix);
