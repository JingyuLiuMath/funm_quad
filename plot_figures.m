function plot_figures(result_list, example_name, save_path)

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
    t = curr_result.param.truncation_length;
    if t ~= inf
        displayname = curr_result.method + " (t = " + string(t) + ")";
    else
        displayname = curr_result.method;
    end
    semilogy(...
        vecnorm(curr_result.f_ex - curr_result.out.appr) / norm(curr_result.f_ex), ...
        '--', ...
        "Marker", marker_list(it), ...
        "Color", color_list(it, :), ...
        "DisplayName",displayname);
    hold on;
end
legend("location", "southeastoutside");
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('rel error');
file_name = example_name + "_rel_err.eps";
saveas(gcf, save_path + file_name, "epsc");

figure();
for it = 1 : m
    curr_result = result_list{it};
    t = curr_result.param.truncation_length;
    if t ~= inf
        displayname = curr_result.method + " (t = " + string(t) + ")";
    else
        displayname = curr_result.method;
    end
    semilogy(...
        curr_result.out.update, ...
        '--', ...
        "Marker", marker_list(it), ...
        "Color", color_list(it, :), ...
        "DisplayName",displayname);
    hold on;
end
legend("location", "southeastoutside");
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('update norm');
file_name = example_name + "_norm_update.eps";
saveas(gcf, save_path + file_name, "epsc");

figure();
for it = 1 : m
    curr_result = result_list{it};
    plot(...
        curr_result.out.num_quadpoints, ...
        '--', ...
        "Marker", marker_list(it), ...
        "Color", color_list(it, :), ...
        "DisplayName",displayname);
    hold on;
end
legend("location", "southeastoutside");
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('num of quad points');
file_name = example_name + "_num_quad.eps";
saveas(gcf, save_path + file_name, "epsc");

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
            "DisplayName",displayname);
        hold on;
    end
end
legend("location", "southeastoutside");
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('subspace dim');
file_name = example_name + "_subspace_dim.eps";
saveas(gcf, save_path + file_name, "epsc");

end