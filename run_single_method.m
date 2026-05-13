function result = run_single_method(A, b, method_name, basic_param, add_param)

rng(1);

switch method_name
    case "FOM"
        param = basic_param;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.arnoldi = "arnoldi";
        param.update = [];

        tic;
        [f, out] = funm_quad_fix(A,b,param);
        t = toc;
    case "sFOM_s"
        param = basic_param;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.sketching_size = add_param.sketching_size;
        param.arnoldi = "sarnoldi";
        param.update = [];

        tic;
        [f,out] = funm_quad_fix(A,b,param);
        t = toc;
    case "adaFOM_t"
        param = basic_param;
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.sketching_size_control = add_param.sketching_size_control;
        param.cond_tol = add_param.cond_tol;
        param.arnoldi = "arnoldi";
        param.update = "last_orth";

        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
    case "sFOM_shm"
        param = basic_param;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.sketching_size = add_param.sketching_size;
        param.arnoldi = "hm-sarnoldi";
        param.update = [];

        tic;
        [f,out] = funm_quad_fix(A,b,param);
        t = toc;
    case "adaFOM_hm_t"
        param = basic_param;
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.sketching_size_control = add_param.sketching_size_control;
        param.cond_tol = add_param.cond_tol;
        param.arnoldi = "arnoldi";
        param.update = "last_hm_orth";

        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
    case "adaFOM_shm_t"
        param = basic_param;
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.sketching_size_control = add_param.sketching_size_control;
        param.cond_tol = add_param.cond_tol;
        param.arnoldi = "arnoldi";
        param.update = "last_hm_sorth";

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

end