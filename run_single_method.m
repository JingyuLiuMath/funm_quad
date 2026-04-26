function result = run_single_method(A, b, method_name, basic_param, add_param)

result = struct();

% fprintf("  function: %s\n", basic_param.function);
% fprintf("  restart_length: %d\n", basic_param.restart_length);
% fprintf("  max_restarts: %d\n", basic_param.max_restarts);

switch method_name
    case "fix_FOM"
        param = basic_param;
        savename = "fix_FOM";

        fprintf("\n\n");
        fprintf("%s\n", savename);

        tic;
        [f, out] = funm_quad_fix(A,b,param);
        t = toc;
    case "fix_sFOM_s"
        param = basic_param;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.sketching_size = add_param.sketching_size;
        param.sarnoldi = 1;
        param.update = "last_sorth";

        savename = "fix_sFOM_s";

        fprintf("\n\n");
        fprintf("%s\n", savename);

        tic;
        [f,out] = funm_quad_fix(A,b,param);
        t = toc;
    case "ada_FOM_t"
        param = basic_param;
        param.max_restarts = add_param.ada_max_restarts;
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.ada_sketching_size_control = add_param.ada_sketching_size_control;
        param.cond_tol = add_param.cond_tol;
        param.sarnoldi = 0;
        param.update = [];

        savename = "ada_FOM_t" + "_" + string(param.truncation_length);

        fprintf("\n\n");
        fprintf("%s\n", savename);

        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
    case "ada_sFOM_t"
        param = basic_param;
        param.max_restarts = add_param.ada_max_restarts;
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.ada_sketching_size_control = add_param.ada_sketching_size_control;
        param.cond_tol = add_param.cond_tol;
        param.sarnoldi = 0;
        param.update = "last_sorth";

        savename = "ada_sFOM_t" + "_" + string(param.truncation_length);

        fprintf("\n\n");
        fprintf("%s\n", savename);

        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
end

result = struct();

result.method = method_name;
result.param = param;
result.f = f;
result.out = out;
result.time = t;
num_it = size(out.appr, 2);
result.num_it = num_it;
result.num_oracle = out.num_oracle;
result.savename = savename;


end