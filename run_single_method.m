function result = run_single_method(A, b, f_ex, method_name, basic_param, add_param)

result = {};

base_add_param = add_param{1};
m = numel(add_param);

% fprintf("\n\n");
% fprintf("%s\n", method_name);
switch method_name
    case "benchmark"
        param = basic_param;
        tic
        [f,out] = funm_quad(A,b,param);
        t = toc;
        result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
    % case "fix FOM-t last-orth"
    %     for k = 1 : m
    %         curr_add_param = add_param{k};
    %         param = basic_param;
    %         param.truncation_length = curr_add_param.truncation_length;
    %         param.max_num_quad_points = curr_add_param.max_num_quad_points;
    %         param.sarnoldi = 0;
    %         param.update = "last_orth";
    %         tic;
    %         [f,out] = funm_quad_fix(A,b,param);
    %         t = toc;
    %         result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
    %     end
    case "fix FOM-t last-sorth"
        for k = 1 : m
            curr_add_param = add_param{k};
            param = basic_param;
            param.truncation_length = curr_add_param.truncation_length;
            param.max_num_quad_points = curr_add_param.max_num_quad_points;
            param.sketching_mat_type = curr_add_param.sketching_mat_type;
            param.sketching_size = curr_add_param.sketching_size;
            param.sarnoldi = 0;
            param.update = "last_sorth";
            tic;
            [f,out] = funm_quad_fix(A,b,param);
            t = toc;
            result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
        end
    % case "fix FOM-t whitening"
    %     for k = 1 : m
    %         curr_add_param = add_param{k};
    %         param = basic_param;
    %         param.truncation_length = curr_add_param.truncation_length;
    %         param.max_num_quad_points = curr_add_param.max_num_quad_points;
    %         param.sketching_mat_type = curr_add_param.sketching_mat_type;
    %         param.sketching_size = curr_add_param.sketching_size;
    %         param.sarnoldi = 0;
    %         param.update = "whitening";
    %         tic;
    %         [f,out] = funm_quad_fix(A,b,param);
    %         t = toc;
    %         result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
    %     end
    case "fix FOM-s"
        param = basic_param;
        param.max_num_quad_points = base_add_param.max_num_quad_points;
        param.sketching_mat_type = base_add_param.sketching_mat_type;
        param.sketching_size = base_add_param.sketching_size;
        param.sarnoldi = 1;
        param.update = "last_sorth";
        tic;
        [f,out] = funm_quad_fix(A,b,param);
        t = toc;
        result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
    case "ada FOM-t last-orth"
        for k = 1 : m
            curr_add_param = add_param{k};
            param = basic_param;
            param.max_restarts = curr_add_param.ada_max_restarts;
            param.truncation_length = curr_add_param.truncation_length;
            param.max_num_quad_points = curr_add_param.max_num_quad_points;
            param.sketching_mat_type = curr_add_param.sketching_mat_type;
            param.ada_sketching_size_control = curr_add_param.ada_sketching_size_control;
            param.cond_tol = curr_add_param.cond_tol;
            param.sarnoldi = 0;
            param.update = "last_orth";
            tic;
            [f,out] = funm_quad_adaptive(A,b,param);
            t = toc;
            result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
        end
    case "ada FOM-t last-sorth"
        for k = 1 : m
            curr_add_param = add_param{k};
            param = basic_param;
            param.max_restarts = curr_add_param.ada_max_restarts;
            param.truncation_length = curr_add_param.truncation_length;
            param.max_num_quad_points = curr_add_param.max_num_quad_points;
            param.sketching_mat_type = curr_add_param.sketching_mat_type;
            param.ada_sketching_size_control = curr_add_param.ada_sketching_size_control;
            param.cond_tol = curr_add_param.cond_tol;
            param.sarnoldi = 0;
            param.update = "last_sorth";
            tic;
            [f,out] = funm_quad_adaptive(A,b,param);
            t = toc;
            result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
        end
    case "ada FOM-t whitening"
        for k = 1 : m
            curr_add_param = add_param{k};
            param = basic_param;
            param.max_restarts = curr_add_param.ada_max_restarts;
            param.truncation_length = curr_add_param.truncation_length;
            param.max_num_quad_points = curr_add_param.max_num_quad_points;
            param.sketching_mat_type = curr_add_param.sketching_mat_type;
            param.ada_sketching_size_control = curr_add_param.ada_sketching_size_control;
            param.cond_tol = curr_add_param.cond_tol;
            param.sarnoldi = 0;
            param.update = "whitening";
            tic;
            [f,out] = funm_quad_adaptive(A,b,param);
            t = toc;
            result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
        end
    case "ada FOM-s last-orth"
        param = basic_param;
        param.max_restarts = base_add_param.ada_max_restarts;
        param.truncation_length = base_add_param.truncation_length;
        param.max_num_quad_points = base_add_param.max_num_quad_points;
        param.sketching_mat_type = base_add_param.sketching_mat_type;
        param.ada_sketching_size_control = base_add_param.ada_sketching_size_control;
        param.cond_tol = base_add_param.cond_tol;
        param.sarnoldi = 1;
        param.update = "last_orth";
        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
        result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
    case "ada FOM-s last-sorth"
        param = basic_param;
        param.max_restarts = base_add_param.ada_max_restarts;
        param.truncation_length = base_add_param.truncation_length;
        param.max_num_quad_points = base_add_param.max_num_quad_points;
        param.sketching_mat_type = base_add_param.sketching_mat_type;
        param.ada_sketching_size_control = base_add_param.ada_sketching_size_control;
        param.cond_tol = base_add_param.cond_tol;
        param.sarnoldi = 1;
        param.update = "last_sorth";
        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
        result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
    case "ada FOM-s whitening"
        param = basic_param;
        param.max_restarts = base_add_param.ada_max_restarts;
        param.truncation_length = base_add_param.truncation_length;
        param.max_num_quad_points = base_add_param.max_num_quad_points;
        param.sketching_mat_type = base_add_param.sketching_mat_type;
        param.ada_sketching_size_control = base_add_param.ada_sketching_size_control;
        param.cond_tol = base_add_param.cond_tol;
        param.sarnoldi = 1;
        param.update = "whitening";
        tic;
        [f,out] = funm_quad_adaptive(A,b,param);
        t = toc;
        result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
    % case "ada FOM-st whitening"
    %     for k = 1 : m
    %         curr_add_param = add_param{k};
    %         param = basic_param;
    %         param.max_restarts = curr_add_param.ada_max_restarts;
    %         param.truncation_length = curr_add_param.truncation_length;
    %         param.max_num_quad_points = curr_add_param.max_num_quad_points;
    %         param.sketching_mat_type = curr_add_param.sketching_mat_type;
    %         param.ada_sketching_size_control = curr_add_param.ada_sketching_size_control;
    %         param.cond_tol = curr_add_param.cond_tol;
    %         param.sarnoldi = 1;
    %         param.update = "whitening";
    %         tic;
    %         [f,out] = funm_quad_adaptive(A,b,param);
    %         t = toc;
    %         result{end + 1} = get_result(f_ex, method_name, param, f, out, t);
    %     end
end

end