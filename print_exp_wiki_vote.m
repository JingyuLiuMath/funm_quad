clear all;
close all;

truncation_length_list = [2, 1, 0];
method_list = [...
    "fix_FOM", ...
    "fix_sFOM_s", ...
    "ada_sFOM_t"
    ];

data_prefix = "./data/wiki_vote_exp/";

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
load('./data/wiki_vote/wiki-Vote.mat');
load('./data/wiki_vote/wiki-Vote-comp.mat');

A = -Problem.A; ee = -ee; N = size(A,1); 
b = ones(N,1);
f_ex = ex_expm;

%% print methods.
caption_name = "Relative error, time, iteration number and oracle number" + ...
    " of the exp function" + ...
    " for the wiki-vote example";
label_name = "exp_wiki_vote";

print_table(data_prefix, f_ex, ...
    method_list, ...
    truncation_length_list, ...
    caption_name, label_name);
