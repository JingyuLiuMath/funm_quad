function plot_figures_thick_restart(data_prefix, f_ex, ...
    method_list, truncation_length_list, ...
    figure_prefix)

results_list_without = [];
results_list_with = [];

num_truncate_len = length(truncation_length_list);
num_method = length(method_list);
for it_method = 1 : num_method
    method_name = method_list(it_method);
    if contains(method_name, "HmArn")
        continue;
    end
    if endsWith(method_name, "_t")
        for it_trunc = 1 : num_truncate_len
            truncation_length = truncation_length_list(it_trunc);
            save_name = method_name + "_" + string(truncation_length);
            load(data_prefix + save_name  + ".mat");
            results_list_without = [results_list_without; result];
            save_name_thick_restart = save_name + "_thick_restart";
            load(data_prefix + save_name_thick_restart  + ".mat");
            results_list_with = [results_list_with; result_with_thick_restart];
        end
    else
        truncation_length = inf;
        save_name = method_name;
        load(data_prefix + save_name  + ".mat");
        results_list_without = [results_list_without; result];
        save_name_thick_restart = save_name + "_thick_restart";
        load(data_prefix + save_name_thick_restart  + ".mat");
        results_list_with = [results_list_with; result_with_thick_restart];
    end
end

plot_result(results_list_without, results_list_with, f_ex, figure_prefix);

end

function plot_result(results_list_without, results_list_with, f_ex, figure_prefix)

num_result = length(results_list_without);

for it = 1 : num_result
    fig = new_plot_figure();
    curr_result_without = results_list_without(it);
    curr_result_with = results_list_with(it);
    save_name = curr_result_without.save_name;
    display_name = replace(save_name, "_", "-");
    fprintf("method: %s\n", display_name);
    err_without = vecnorm(f_ex - curr_result_without.out.appr) / norm(f_ex);
    err_with = vecnorm(f_ex - curr_result_with.out.appr) / norm(f_ex);

    semilogy(cumsum(curr_result_without.num_oracle), err_without, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Marker", "+", ...
        'Color', "g", ...
        "MarkerSize", 10);
    total_oracle_without = sum(curr_result_without.num_oracle);
    fprintf("number of oracle calls without thick restart: %d\n", total_oracle_without);
    fprintf("final relative error without thick restart: %.2e\n", err_without(end));
    hold on;
    semilogy(cumsum(curr_result_with.num_oracle), err_with, ...
        "LineWidth", 2, ...
        "LineStyle", "--", ...
        "Marker", "x", ...
        'Color', "m", ...
        "MarkerSize", 10);
    total_oracle_with = sum(curr_result_with.num_oracle);
    fprintf("number of oracle calls with thick restart: %d\n", total_oracle_with);
    fprintf("final relative error with thick restart: %.2e\n", err_with(end));
    reduced_ratio = (total_oracle_without - total_oracle_with) / total_oracle_without;
    fprintf("reduced ratio: %.2f%%\n", reduced_ratio * 100);

    % if length(err_without) <= 10 &&  length(err_with) <= 10
    %     for k = 2 : length(err_without)
    %         text(k + 0.1, 2 * err_without(k), num2str(curr_result_without.out.num_quadpoints(k)), ...
    %             'Color', [0 1 0], 'FontSize', 16, 'Rotation', 45);
    %     end
    % 
    %     for k = 2 : length(err_with)
    %         text(k + 0.1, 2 * err_with(k), num2str(curr_result_with.out.num_quadpoints(k)), ...
    %             'Color', [1 0 1], 'FontSize', 16, 'Rotation', 45);
    %     end
    % end
    display_name_without = "without thick restart (" + sprintf('%.1f', curr_result_without.time) + ")";
    display_name_with = "with thick restart (" + sprintf('%.1f', curr_result_with.time) + ")";
    legend(display_name_without, display_name_with);
    xlabel('number of matrix-vector products');
    ylabel('rel error');
    title(display_name);
    set(gca, 'FontSize', 40);
    hold off;
    save_plot_figure(fig, figure_prefix + "_" + display_name + "_thick_restart_rel_err.pdf", "pdf")
    save_plot_figure(fig, figure_prefix + "_" + display_name + "_thick_restart_rel_err.png", "png");
end

end

function fig = new_plot_figure()

figure_size_pixels = [1892, 1028];
figure_resolution = 100;
fig = figure(...
    "Units", "pixels", ...
    "Position", [100, 100, figure_size_pixels]);
paper_size_inches = figure_size_pixels / figure_resolution;
set(fig, ...
    "PaperUnits", "inches", ...
    "PaperPositionMode", "manual", ...
    "PaperPosition", [0, 0, paper_size_inches], ...
    "PaperSize", paper_size_inches);

end

function save_plot_figure(fig, file_name, file_format)

drawnow;
switch file_format
    case "epsc"
        print(fig, char(file_name), "-depsc", "-painters");
    case "pdf"
        print(fig, char(file_name), "-dpdf", "-painters");
    case "png"
        print(fig, char(file_name), "-dpng", "-r100");
    otherwise
        saveas(fig, file_name, file_format);
end

end
