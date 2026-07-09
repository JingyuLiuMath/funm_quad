function run_methods_thick_restart(data_prefix, ...
    method_list, ...
    truncation_length_list, ...
    basic_param, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    sketching_size_control, cond_tol, ...
    A, b, ...
    number_thick)

num_truncate_len = length(truncation_length_list);
num_method = length(method_list);
for it_method = 1 : num_method
    method_name = method_list(it_method);
    if contains(method_name, "GMRES")
        continue;
    end
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

            thick_param = basic_param;
            thick_param.thick = @thick_quad;
            thick_param.number_thick = number_thick;

            result_with_thick_restart = run_single_method(A, b, method_name, thick_param, add_param);
            save_name_with = save_name;
            result_with_thick_restart.save_name = save_name_with;

            save(data_prefix + save_name + "_thick_restart" + ".mat", "result_with_thick_restart");
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

        thick_param = basic_param;
        thick_param.thick = @thick_quad;
        thick_param.number_thick = number_thick;

        result_with_thick_restart = run_single_method(A, b, method_name, thick_param, add_param);
        save_name_with = save_name + " (with thick restarting)";
        result_with_thick_restart.save_name = save_name_with;

        save(data_prefix + save_name + "_thick_restart" + ".mat", "result_with_thick_restart");
    end
end

end