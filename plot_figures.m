function plot_figures(result_list, save_path, truncation_length)

m = length(result_list);

color_list = get_colors(m);
marker_list = get_markers(m);

max_iter = 0;
for it = 1 : m
    max_iter = max(max_iter, result_list{it}.num_it);
end

figure();
for it = 1 : m
    curr_result = result_list{it};
    semilogy(...
        vecnorm(curr_result.f_ex - curr_result.out.appr) / norm(curr_result.f_ex), ...
        '--', ...
        "Marker", marker_list(it), ...
        "Color", color_list(it, :), ...
        "DisplayName", curr_result.method);
    hold on;
end
legend;
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('rel error');
file_name = "wikivote_rel_err_trunc_" + string(truncation_length);
saveas(gcf, save_path + file_name + ".eps", "epsc");

figure();
for it = 1 : m
    curr_result = result_list{it};
    semilogy(...
        curr_result.out.update, ...
        '--', ...
        "Marker", marker_list(it), ...
        "Color", color_list(it, :), ...
        "DisplayName", curr_result.method);
    hold on;
end
legend;
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('update norm');
file_name = "wikivote_norm_update_trunc_" + string(truncation_length);
saveas(gcf, save_path + file_name + ".eps", "epsc");
saveas(gcf, "./figure/wikivote/" + file_name + ".eps", "epsc");

figure();
for it = 1 : m
    curr_result = result_list{it};
    plot(...
        curr_result.out.num_quadpoints, ...
        '--', ...
        "Marker", marker_list(it), ...
        "Color", color_list(it, :), ...
        "DisplayName", curr_result.method);
    hold on;
end
legend;
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('num of quad points');
file_name = "wikivote_num_quad_trunc_" + string(truncation_length);
saveas(gcf, save_path + file_name + ".eps", "epsc");

figure();
for it = 1 : m
    curr_result = result_list{it};
    prefix_char = char(curr_result.method);
    ada_flag = strcmp(prefix_char(1:3), 'ada');
    if ada_flag
        plot(curr_result.out.dim, ...
            '--', ...
            "Marker", marker_list(it), ...
            "Color", color_list(it, :), ...
            "DisplayName", curr_result.method);
        hold on;
    end
end
legend;
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('subspace dim');
file_name = "wikivote_subspace_dim_trunc_" + string(truncation_length);
saveas(gcf, save_path + file_name + ".eps", "epsc");

end