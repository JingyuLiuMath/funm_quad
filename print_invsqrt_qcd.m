clear;
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
fprintf("file_prefix: %s\n", data_prefix);

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
fprintf("\n");
fprintf("\\begin{table}[tbhp]\n")
fprintf("\\centering\n")
fprintf("\\begin{tabular}{ccccc}\n")
fprintf("\\toprule\n")
fprintf("Method & error & time & iter number & oracle number \\\\ \n")
num_truncate_len = length(truncation_length_list);
num_method = length(method_list);
for it_method = 1 : num_method
    method_name = method_list(it_method);
    if endsWith(method_name, "_t")
        for it_trunc = 1 : num_truncate_len
            truncation_length = truncation_length_list(it_trunc);
            save_name = method_name + "_" + string(truncation_length);
            load(data_prefix + save_name  + ".mat");
            print_result(save_name, result, f_ex);
        end
    else
        truncation_length = inf;
        save_name = method_name;
        load(data_prefix + save_name  + ".mat");
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
fprintf("& %d ", result.num_it);
fprintf("& %d ", result.num_oracle);
fprintf("\\\\ \n");

end
