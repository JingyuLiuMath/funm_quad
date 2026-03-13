function result_list = run_methods(A, b, f_ex, method_list, basic_param, add_param)

m = length(method_list);
result_list = cell(1, m);

for it = 1 : m
    curr_method = method_list(it);
    result_list{it} = run_single_method(A, b, f_ex, curr_method, basic_param, add_param);
end

end