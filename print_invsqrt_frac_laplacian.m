clear all;
close all;

data_prefix = "./data/frac_laplacian_invsqrt/";

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

%% Load matrix.
load('./data/frac_laplacian/gnutella_comp.mat');
f_ex = ex_gnutella;

%% print methods.
caption_name = "Relative error, time, iteration number and oracle number" + ...
    " of the invsqrt function" + ...
    " for the fractional Laplacian example";
label_name = "invsqrt_frac_laplacian";

print_table(data_prefix, f_ex, ...
    method_list, ...
    truncation_length_list, ...
    caption_name, label_name);
