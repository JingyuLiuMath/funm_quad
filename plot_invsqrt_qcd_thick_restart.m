clear all;
close all;

truncation_length_list = [2, 1, 0];
method_list = [...
    "fix_FOM", ...
    "fix_sFOM_s", ...
    "ada_sFOM_t"
    ];

data_prefix = "./data/qcd_invsqrt_thick_restart/";
figure_prefix = "./figure/qcd_invsqrt_thick_restart/qcd_invsqrt_thick_restart";

results_list_without = [];
results_list_with = [];

fprintf("truncation_length: ");
for it = 1 : length(truncation_length_list)
    fprintf("%d ", truncation_length_list(it));
end
fprintf("\n");
fprintf("file_prefix: %s\n", data_prefix);

%% Load matrix.
load("./data/qcd/qcd_matrix_nonhermitian.mat");
load("./data/qcd/qcd_nonhermitian_exact.mat");
c = zeros(size(Q,1),1); c(1) = 1;
b = Q*c; normv = norm(b);
f_ex = qcd_nonhermitian_exact / normv;

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
            results_list_without = [results_list_without; result_without_thick_restart];
            results_list_with = [results_list_with; result_with_thick_restart];
        end
    else
        truncation_length = inf;
        save_name = method_name;
        load(data_prefix + save_name  + ".mat");
        results_list_without = [results_list_without; result_without_thick_restart];
        results_list_with = [results_list_with; result_with_thick_restart];
    end
end

plot_result(results_list_without, results_list_with, f_ex, figure_prefix);

function plot_result(results_list_without, results_list_with, f_ex, figure_prefix)

num_result = length(results_list_without);

color_list = get_colors(num_result);
marker_list = get_markers(num_result);

for it = 1 : num_result
    figure();
    curr_result_without = results_list_without(it);
    curr_result_with = results_list_with(it);
    save_name = curr_result_without.save_name;
    display_name = replace(save_name, "_", "-");

    err_without = vecnorm(f_ex - curr_result_without.out.appr) / norm(f_ex);
    err_with = vecnorm(f_ex - curr_result_with.out.appr) / norm(f_ex);

    semilogy(err_without, 'g--+');
    hold on;
    semilogy(err_with, 'm--+');

    if length(err_without) <= 10 &&  length(err_with) <= 10
        for k = 2 : length(err_without)
            text(k + 0.1, 2 * err_without(k), num2str(curr_result_without.out.num_quadpoints(k)), ...
                'Color', [0 1 0], 'FontSize', 16, 'Rotation', 45);
        end

        for k = 2 : length(err_with)
            text(k + 0.1, 2 * err_with(k), num2str(curr_result_with.out.num_quadpoints(k)), ...
                'Color', [1 0 1], 'FontSize', 16, 'Rotation', 45);
        end
    end
    display_name_without = "without thick restart (" + num2str(curr_result_without.time) + ")";
    display_name_with = "with thick restart (" + num2str(curr_result_with.time) + ")";
    legend(display_name_without, display_name_with);
    xlabel('cycle');
    ylabel('rel err');
    title(display_name);
    hold off;
    saveas(gcf, figure_prefix + "_" + save_name + "_rel_err.eps", "epsc")
end

end
