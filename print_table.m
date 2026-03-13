function print_table(result_list)

m = length(result_list);

for it = 1 : m
    curr_result = result_list{it};
    fprintf("%s & %d & %.4e & %.4e \\\\ \n", ...
        curr_result.method, curr_result.num_it, curr_result.rel_err, curr_result.t);
end

end