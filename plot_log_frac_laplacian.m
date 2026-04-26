clear;
close all;

truncation_length_list = [2, 1, 0];
method_list = [...
    "fix_FOM", ...
    "fix_sFOM_s", ...
    "ada_sFOM_t"
    ];
data_prefix = "./data/frac_laplacian_log/";
figure_prefix = "./figure/frac_laplacian_log/frac_laplacian_log";
results_list = [];


fprintf("truncation_length: ");
for it = 1 : length(truncation_length_list)
    fprintf("%d ", truncation_length_list(it));
end
fprintf("\n");
fprintf("file_prefix: %s\n", data_prefix);

%% Load matrix.
load('./data/frac_laplacian/gnutella_comp.mat');
A = L;
b = A * b;

addpath('funm_quad')
exact_param.function = 'log';
exact_param.restart_length = 300; 
exact_param.max_restarts = 15;
exact_param.tol = 1e-15;
exact_param.transformation_parameter = 1;
exact_param.hermitian = 0;
exact_param.V_full = 0;
exact_param.H_full = 0;
exact_param.exact = [];
exact_param.stopping_accuracy = 1e-15;
exact_param.inner_product = @(a,b) b'*a;
exact_param.thick = [];
exact_param.min_decay = .95;
exact_param.waitbar = 0;
exact_param.reorth_number = 0;
exact_param.truncation_length = inf;
exact_param.verbose = 1;

[f_ex, ~] = funm_quad(A, b, exact_param);

%% print methods.
num_truncate_len = length(truncation_length_list);
num_method = length(method_list);
for it_method = 1 : num_method
    method_name = method_list(it_method);
    if endsWith(method_name, "_t")
        for it_trunc = 1 : num_truncate_len
            truncation_length = truncation_length_list(it_trunc);
            save_name = method_name + "_" + string(truncation_length);
            load(data_prefix + save_name  + ".mat");
            results_list = [results_list; result];
        end
    else
        truncation_length = inf;
        save_name = method_name;
        load(data_prefix + save_name  + ".mat");
        results_list = [results_list; result];
    end
end

plot_result(results_list, f_ex, figure_prefix);

function plot_result(results_list, f_ex, figure_prefix)

num_result = length(results_list);

color_list = get_colors(num_result);
marker_list = get_markers(num_result);

figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");

    loglog(vecnorm(curr_result.out.appr - f_ex) / norm(f_ex), ...
        '--', ...
        "Color", color_list(it, :), ...
        "Marker", marker_list(it, :), ...
        "DisplayName", display_name);
    hold on;
end
legend("location", "southeastoutside");
xlabel('cycle');
ylabel('rel error');
file_name = figure_prefix + "_rel_err.eps";
saveas(gcf, file_name, "epsc");

figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");

    loglog(curr_result.out.update, ...
        '--', ...
        "Color", color_list(it, :), ...
        "Marker", marker_list(it, :), ...
        "DisplayName", display_name);
    hold on;
end
legend("location", "southeastoutside");
xlabel('cycle');
ylabel('update norm');
file_name = figure_prefix + "_norm_update.eps";
saveas(gcf, file_name, "epsc");

figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");

    plot(curr_result.out.num_quadpoints, ...
        '--', ...
        "Color", color_list(it, :), ...
        "Marker", marker_list(it, :), ...
        "DisplayName", display_name);
    hold on;
end
legend("location", "southeastoutside");
xlabel('cycle');
ylabel('num of quad points');
file_name = figure_prefix + "_num_quad.eps";
saveas(gcf, file_name, "epsc");

figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");
    
    prefix_char = char(curr_result.method);
    ada_flag = strcmp(prefix_char(1:3), 'ada');
    if ada_flag
        plot(curr_result.out.dim, ...
            '--', ...
            "Color", color_list(it, :), ...
            "Marker", marker_list(it, :), ...
            "DisplayName", display_name);
        hold on;
    end
end
legend("location", "southeastoutside");
xlabel('cycle');
ylabel('subspace dim');
file_name = figure_prefix + "_subspace_dim.eps";
saveas(gcf, file_name, "epsc");

end
