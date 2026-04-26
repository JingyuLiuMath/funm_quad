function add_param = construct_ada_param(...
    truncation_length, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    ada_max_restarts, ada_sketching_size_control, cond_tol)

add_param = struct();
add_param.truncation_length = truncation_length;
add_param.max_num_quad_points = max_num_quad_points;
add_param.sketching_mat_type = sketching_mat_type;
add_param.sketching_size = sketching_size;
add_param.ada_max_restarts = ada_max_restarts;
add_param.ada_sketching_size_control = ada_sketching_size_control;
add_param.cond_tol = cond_tol;

end