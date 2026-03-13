function print_table(result_list)

m = length(result_list);

for it = 1 : m
    curr_result = result_list{it};
    t = curr_result.param.truncation_length;
    if t ~= inf
        fprintf("%s (t = %d) & %d & %.4e & %.4e \\\\ \n", ...
            curr_result.method, t, curr_result.num_it, curr_result.rel_err, curr_result.time);
    else
        fprintf("%s & %d & %.4e & %.4e \\\\ \n", ...
            curr_result.method, curr_result.num_it, curr_result.rel_err, curr_result.time);
    end
end

end