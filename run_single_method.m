function result = run_single_method(A, b, f_ex, method_name, basic_param, add_param)

param = basic_param;

fprintf("\n\n");
fprintf("%s\n", method_name);
switch method_name
    case "benchmark"
        tic
        [f,out] = funm_quad(A,b,param);
        t = toc;
    case "fix FOM-t last-orth"
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sarnoldi = 0;
        param.update = "last_orth";
        tic;
        [f,out] = funm_quad_fix(A,b,param);
        t = toc;
    case "fix FOM-t last-sorth"
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.sketching_size = add_param.sketching_size;
        param.sarnoldi = 0;
        param.update = "last_sorth";
        tic;
        [f,out] = funm_quad_fix(A,b,param);
        t = toc;
    case "fix FOM-t whitening"
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.sketching_size = add_param.sketching_size;
        param.sarnoldi = 0;
        param.update = "whitening";
        tic;
        [f,out] = funm_quad_fix(A,b,param);
        t = toc;
    case "fix FOM-s"
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.sketching_size = add_param.sketching_size;
        param.sarnoldi = 1;
        param.update = [];
        tic;
        [f,out] = funm_quad_fix(A,b,param);
        t = toc;
    case "ada FOM-t last-orth"
        param.max_restarts = add_param.ada_max_restarts;
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.ada_sketching_size_control = add_param.ada_sketching_size_control;
        param.cond_tol = add_param.cond_tol;
        param.sarnoldi = 0;
        param.update = "last_orth";
        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
    case "ada FOM-t last-sorth"
        param.max_restarts = add_param.ada_max_restarts;
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.ada_sketching_size_control = add_param.ada_sketching_size_control;
        param.cond_tol = add_param.cond_tol;
        param.sarnoldi = 0;
        param.update = "last_sorth";
        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
    case "ada FOM-t whitening"
        param.max_restarts = add_param.ada_max_restarts;
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.ada_sketching_size_control = add_param.ada_sketching_size_control;
        param.cond_tol = add_param.cond_tol;
        param.sarnoldi = 0;
        param.update = "whitening";
        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
    case "ada FOM-s whitening"
        param.max_restarts = add_param.ada_max_restarts;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.ada_sketching_size_control = add_param.ada_sketching_size_control;
        param.cond_tol = add_param.cond_tol;
        param.sarnoldi = 1;
        param.update = "whitening";
        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
    case "ada FOM-st whitening"
        param.max_restarts = add_param.ada_max_restarts;
        param.truncation_length = add_param.truncation_length;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.ada_sketching_size_control = add_param.ada_sketching_size_control;
        param.cond_tol = add_param.cond_tol;
        param.sarnoldi = 1;
        param.update = "whitening";
        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
end

num_it = size(out.appr, 2);
rel_err = norm(f_ex - f) / norm(f_ex);

result = struct();
result.method = method_name;
result.f = f;
result.out = out;
result.num_it = num_it;
result.rel_err = rel_err;
result.t = t;
result.f_ex = f_ex;
result.param = param;

fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t)
fprintf("\n\n");
