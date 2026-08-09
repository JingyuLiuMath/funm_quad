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
            result.save_name = method_name + " (t = " + string(truncation_length) + ")";
            results_list = [results_list; result];
        end
    else
        truncation_length = inf;
        save_name = method_name;
        load(data_prefix + save_name  + ".mat");
        result.save_name = method_name;
        results_list = [results_list; result];
    end
end

plot_result(results_list, f_ex, figure_prefix);

end

function plot_result(results_list, f_ex, figure_prefix)

num_result = length(results_list);

color_list = get_colors(num_result);
marker_list = get_markers(num_result);

fig = new_plot_figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");

    semilogy(cumsum(curr_result.num_oracle), vecnorm(curr_result.out.appr - f_ex) / norm(f_ex), ...
        '--', ...
        "Color", color_list(it, :), ...
        "Marker", marker_list(it, :), ...
        "LineWidth", 4, ...
        "MarkerSize", 20, ...
        "DisplayName", display_name);
    hold on;
end
legend("location", "southeastoutside");
xlabel('number of matrix-vector products');
ylabel('relative error');
set(gca, 'FontSize', 40);
file_name = figure_prefix + "_rel_err.pdf";
save_plot_figure(fig, file_name, "pdf");
file_name = figure_prefix + "_rel_err.png";
save_plot_figure(fig, file_name, "png");

fig = new_plot_figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");

    semilogy(curr_result.out.update, ...
        '--', ...
        "Color", color_list(it, :), ...
        "Marker", marker_list(it, :), ...
        "LineWidth", 4, ...
        "MarkerSize", 20, ...
        "DisplayName", display_name);
    hold on;
end
legend("location", "northeast");
xlabel('iteration');
ylabel('update norm');
set(gca, 'FontSize', 40);
file_name = figure_prefix + "_norm_update.pdf";
save_plot_figure(fig, file_name, "pdf");
file_name = figure_prefix + "_norm_update.png";
save_plot_figure(fig, file_name, "png");

fig = new_plot_figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");

    plot(curr_result.out.num_quadpoints, ...
        '--', ...
        "Color", color_list(it, :), ...
        "Marker", marker_list(it, :), ...
        "LineWidth", 4, ...
        "MarkerSize", 20, ...
        "DisplayName", display_name);
    hold on;
end
legend("location", "northeast");
xlabel('iteration');
ylabel('num of quad points');
set(gca, 'FontSize', 40);
file_name = figure_prefix + "_num_quad.pdf";
save_plot_figure(fig, file_name, "pdf");
file_name = figure_prefix + "_num_quad.png";
save_plot_figure(fig, file_name, "png");

fig = new_plot_figure();
for it = 1 : num_result
    curr_result = results_list(it);
    display_name = replace(curr_result.save_name, "_", "-");
    fprintf("method: %s, subspace dim: %d\n", display_name, round(mean(curr_result.out.dim)));
    plot(curr_result.out.dim, ...
            '--', ...
            "Color", color_list(it, :), ...
            "Marker", marker_list(it, :), ...
            "LineWidth", 4, ...
            "MarkerSize", 20, ...
            "DisplayName", display_name);
        hold on;
end
legend("location", "northeast");
xlabel('iteration');
ylabel('subspace dim');
set(gca, 'FontSize', 40);
file_name = figure_prefix + "_subspace_dim.pdf";
save_plot_figure(fig, file_name, "pdf");
file_name = figure_prefix + "_subspace_dim.png";
save_plot_figure(fig, file_name, "png");

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
