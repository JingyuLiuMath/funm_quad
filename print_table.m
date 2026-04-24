function print_table(result_list)

m = length(result_list);

fprintf("\\begin{table}[tbhp]\n")
fprintf("\\centering\n")
fprintf("\\begin{tabular}{cccc}\n")
fprintf("\\toprule\n")
fprintf("Method & \\(n_{\\iter}\\) & \\(e\\) & \\(t_{\\iter}\\) \\\\ \n")
for it = 1 : m
    fprintf("\\midrule\n")
    curr_result = result_list{it};
    t = curr_result.param.truncation_length;
    if t ~= inf
        fprintf("%s (t = %d) & %d & %.1e & %.1e \\\\ \n", ...
            curr_result.method, t, curr_result.num_it, curr_result.rel_err, curr_result.time);
    else
        fprintf("%s & %d & %.1e & %.1e \\\\ \n", ...
            curr_result.method, curr_result.num_it, curr_result.rel_err, curr_result.time);
    end
end
fprintf("\\bottomrule\n")
fprintf("\\end{tabular}\n")

fprintf("\\end{table}\n")
fprintf("\n\n")

end