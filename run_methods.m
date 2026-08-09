function run_methods(data_prefix, ...
    method_list, ...
    truncation_length_list, ...
    basic_param, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    sketching_size_control, cond_tol, ...
    A, b)

num_truncate_len = length(truncation_length_list);
num_method = length(method_list);
for it_method = 1 : num_method
    method_name = method_list(it_method);
    if endsWith(method_name, "_t")
        for it_trunc = 1 : num_truncate_len
            truncation_length = truncation_length_list(it_trunc);
            
            save_name = method_name + "_" + string(truncation_length);
            fprintf("\n\n");
            fprintf("%s\n", save_name);

            add_param = construct_ada_param(...
                truncation_length, ...
                max_num_quad_points, ...
                sketching_mat_type, sketching_size, ...
                sketching_size_control, cond_tol);

            result = run_single_method(A, b, method_name, basic_param, add_param);
            result.save_name = save_name;
            
            save(data_prefix + save_name  + ".mat", "result");
        end
    else
        truncation_length = inf;

        save_name = method_name;
        fprintf("\n\n");
        fprintf("%s\n", save_name);
        
        add_param = construct_ada_param(...
            truncation_length, ...
            max_num_quad_points, ...
            sketching_mat_type, sketching_size, ...
            sketching_size_control, cond_tol);

        result = run_single_method(A, b, method_name, basic_param, add_param);
        result.save_name = save_name;
        
        save(data_prefix + save_name  + ".mat", "result");
    end
end

end