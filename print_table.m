function print_table(data_prefix, f_ex, ...
    method_list, ...
    truncation_length_list, ...
    caption_name, label_name)

fprintf("\n");
fprintf("\\begin{table}[tbhp]\n")
fprintf("\\centering\n")
fprintf("\\caption{%s}\n", caption_name);
fprintf("\\label{tab:%s}\n", label_name);
fprintf("\\begin{tabular}{cccccc}\n")
fprintf("\\toprule\n")
fprintf("method & error & time (s) & \\(n_{\\iter}\\) & \\(n_{\\matvec}\\) & \\(n_{\\vecs}\\)\\\\ \n")
num_truncate_len = length(truncation_length_list);
num_method = length(method_list);
for it_method = 1 : num_method
    method_name = method_list(it_method);
    if endsWith(method_name, "_t")
        for it_trunc = 1 : num_truncate_len
            truncation_length = truncation_length_list(it_trunc);
            save_name = method_name + "_" + string(truncation_length);
            print_name = method_name + " (t = " + string(truncation_length) + ")";
            load(data_prefix + save_name  + ".mat");
            print_result(print_name, result, f_ex);
        end
    else
        truncation_length = inf;
        save_name = method_name;
        print_name = save_name;
        load(data_prefix + save_name  + ".mat");
        print_result(print_name, result, f_ex);
    end
end
fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\end{table}\n");
fprintf("\n");

end

function print_result(save_name, result, f_ex)

rel_err = norm(f_ex - result.f) / norm(f_ex);
print_name = replace(save_name, "_", "-");
fprintf("\\midrule\n")
fprintf("%s ", print_name );
fprintf("& %.1e ", rel_err);
fprintf("& %.1e ", result.time);
fprintf("& %d ", result.num_it);
fprintf("& %d ", sum(result.num_oracle));
fprintf("& %d ", max(result.out.dim));
fprintf("\\\\ \n");

end