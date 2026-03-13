function result = run_single_method(A, b, f_ex, method_name, basic_param, add_param)

result = {};

m = length(add_param.truncation_length);

fprintf("\n\n");
fprintf("%s\n", method_name);
switch method_name
    case "benchmark"
        param = basic_param;
        tic
        [f,out] = funm_quad(A,b,param);
        t = toc;
        result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
    case "fix FOM-t last-orth"
        for k = 1 : m
            param = basic_param;
            param.truncation_length = add_param.truncation_length(k);
            param.max_num_quad_points = add_param.max_num_quad_points;
            param.sarnoldi = 0;
            param.update = "last_orth";
            tic;
            [f,out] = funm_quad_fix(A,b,param);
            t = toc;
            result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
        end
    case "fix FOM-t last-sorth"
        for k = 1 : m
            param = basic_param;
            param.truncation_length = add_param.truncation_length(k);
            param.max_num_quad_points = add_param.max_num_quad_points;
            param.sketching_mat_type = add_param.sketching_mat_type;
            param.sketching_size = add_param.sketching_size;
            param.sarnoldi = 0;
            param.update = "last_sorth";
            tic;
            [f,out] = funm_quad_fix(A,b,param);
            t = toc;
            result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
        end
    case "fix FOM-t whitening"
        for k = 1 : m
            param = basic_param;
            param.truncation_length = add_param.truncation_length(k);
            param.max_num_quad_points = add_param.max_num_quad_points;
            param.sketching_mat_type = add_param.sketching_mat_type;
            param.sketching_size = add_param.sketching_size;
            param.sarnoldi = 0;
            param.update = "whitening";
            tic;
            [f,out] = funm_quad_fix(A,b,param);
            t = toc;
            result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
        end
    case "fix FOM-s"
        param = basic_param;
        param.max_num_quad_points = add_param.max_num_quad_points;
        param.sketching_mat_type = add_param.sketching_mat_type;
        param.sketching_size = add_param.sketching_size;
        param.sarnoldi = 1;
        param.update = [];
        tic;
        [f,out] = funm_quad_fix(A,b,param);
        t = toc;
        result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
    case "ada FOM-t last-orth"
        for k = 1 : m
            param = basic_param;
            param.max_restarts = add_param.ada_max_restarts;
            param.truncation_length = add_param.truncation_length(k);
            param.max_num_quad_points = add_param.max_num_quad_points;
            param.sketching_mat_type = add_param.sketching_mat_type;
            param.ada_sketching_size_control = add_param.ada_sketching_size_control;
            param.cond_tol = add_param.cond_tol;
            param.sarnoldi = 0;
            param.update = "last_orth";
            tic;
            [f,out] = funm_quad_adaptive(A,b,param);
            t = toc;
            result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
        end
    case "ada FOM-t last-sorth"
        for k = 1 : m
            param = basic_param;
            param.max_restarts = add_param.ada_max_restarts;
            param.truncation_length = add_param.truncation_length(k);
            param.max_num_quad_points = add_param.max_num_quad_points;
            param.sketching_mat_type = add_param.sketching_mat_type;
            param.ada_sketching_size_control = add_param.ada_sketching_size_control;
            param.cond_tol = add_param.cond_tol;
            param.sarnoldi = 0;
            param.update = "last_sorth";
            tic;
            [f,out] = funm_quad_adaptive(A,b,param);
            t = toc;
            result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
        end
    case "ada FOM-t whitening"
        for k = 1 : m
            param = basic_param;
            param.max_restarts = add_param.ada_max_restarts;
            param.truncation_length = add_param.truncation_length(k);
            param.max_num_quad_points = add_param.max_num_quad_points;
            param.sketching_mat_type = add_param.sketching_mat_type;
            param.ada_sketching_size_control = add_param.ada_sketching_size_control;
            param.cond_tol = add_param.cond_tol;
            param.sarnoldi = 0;
            param.update = "whitening";
            tic;
            [f,out] = funm_quad_adaptive(A,b,param);
            t = toc;
            result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
        end
    case "ada FOM-s whitening"
        param = basic_param;
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
        result = get_result(f_ex, method_name, param, f, out, t);
    case "ada FOM-st whitening"
        for k = 1 : m
            param = basic_param;
            param.max_restarts = add_param.ada_max_restarts;
            param.truncation_length = add_param.truncation_length(k);
            param.max_num_quad_points = add_param.max_num_quad_points;
            param.sketching_mat_type = add_param.sketching_mat_type;
            param.ada_sketching_size_control = add_param.ada_sketching_size_control;
            param.cond_tol = add_param.cond_tol;
            param.sarnoldi = 1;
            param.update = "whitening";
            tic;
            [f,out] = funm_quad_adaptive(A,b,param);
            t = toc;
            result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
        end
end