function plot_figures(data_prefix, f_ex, ...
    method_list, truncation_length_list, ...
    figure_prefix)

results_list = [];
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

end

function plot_result(results_list, f_ex, figure_prefix)

num_result = length(results_list);

color_list = get_colors(num_result);
marker_list = get_markers(num_result);

figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");

    semilogy(cumsum(curr_result.num_oracle), vecnorm(curr_result.out.appr - f_ex) / norm(f_ex), ...
        '--', ...
        "Color", color_list(it, :), ...
        "Marker", marker_list(it, :), ...
        "LineWidth", 2, ...
        "DisplayName", display_name);
    hold on;
end
legend("location", "southeastoutside");
xlabel('number of matrix-vector products');
ylabel('relative error');
set(gca, 'FontSize', 18);
file_name = figure_prefix + "_rel_err.eps";
saveas(gcf, file_name, "epsc");
file_name = figure_prefix + "_rel_err.png";
saveas(gcf, file_name, "png");

figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");

    semilogy(curr_result.out.update, ...
        '--', ...
        "Color", color_list(it, :), ...
        "Marker", marker_list(it, :), ...
        "LineWidth", 2, ...
        "DisplayName", display_name);
    hold on;
end
legend("location", "southeastoutside");
xlabel('iteration');
ylabel('update norm');
set(gca, 'FontSize', 18);
file_name = figure_prefix + "_norm_update.eps";
saveas(gcf, file_name, "epsc");
file_name = figure_prefix + "_norm_update.png";
saveas(gcf, file_name, "png");

figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");

    plot(curr_result.out.num_quadpoints, ...
        '--', ...
        "Color", color_list(it, :), ...
        "Marker", marker_list(it, :), ...
        "LineWidth", 2, ...
        "DisplayName", display_name);
    hold on;
end
legend("location", "southeastoutside");
xlabel('iteration');
ylabel('num of quad points');
set(gca, 'FontSize', 18);
file_name = figure_prefix + "_num_quad.eps";
saveas(gcf, file_name, "epsc");
file_name = figure_prefix + "_num_quad.png";
saveas(gcf, file_name, "png");

figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");
    
    plot(curr_result.out.dim, ...
            '--', ...
            "Color", color_list(it, :), ...
            "Marker", marker_list(it, :), ...
            "LineWidth", 2, ...
            "DisplayName", display_name);
        hold on;
end
legend("location", "southeastoutside");
xlabel('iteration');
ylabel('subspace dim');
set(gca, 'FontSize', 18);
file_name = figure_prefix + "_subspace_dim.eps";
saveas(gcf, file_name, "epsc");
file_name = figure_prefix + "_subspace_dim.png";
saveas(gcf, file_name, "png");

end