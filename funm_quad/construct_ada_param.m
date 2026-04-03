function add_param = construct_ada_param(...
    truncation_length, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    ada_max_restarts, ada_sketching_size_control, cond_tol)

truncation_values = truncation_length(:).';
if isempty(truncation_values)
    truncation_values = inf;
end

add_param = cell(1, numel(truncation_values));
for k = 1 : numel(truncation_values)
    add_param{k} = make_single_add_param(truncation_values(k), ...
        max_num_quad_points, sketching_mat_type, sketching_size, ...
        ada_max_restarts, ada_sketching_size_control, cond_tol);
end

end

function item = make_single_add_param(...
    truncation_length, max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    ada_max_restarts, ada_sketching_size_control, cond_tol)

item = struct();
item.truncation_length = truncation_length;
item.max_num_quad_points = max_num_quad_points;
item.sketching_mat_type = sketching_mat_type;
item.sketching_size = sketching_size;
item.ada_max_restarts = ada_max_restarts;
item.ada_sketching_size_control = ada_sketching_size_control;
item.cond_tol = cond_tol;

end