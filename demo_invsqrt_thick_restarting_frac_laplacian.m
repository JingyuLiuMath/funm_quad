clear;
close all;
rng(2026);
maxNumCompThreads(1);
warning off;

quad_tol = 1e-7;
stop_tol = 1e-8;

m = 300;
max_restarts = 15;
ada_max_restarts = 200;
truncation_length = 2;
max_num_quad_points = 8192;

sketching_mat_type = "sparse sign";
sketching_size = 2 * m;
ada_sketching_size_control = 2;
cond_tol = 1e4;

check_flag = 1;

method_list = [...
    "benchmark", ...
    "fix FOM-s", ...
    "ada FOM-t last-sorth"
    ];
method_name = "ada FOM-t last-sorth";

fprintf("quad_tol: %.1e\n", quad_tol);
fprintf("stop_tol: %.1e\n", stop_tol);

fprintf("m: %d\n", m);
fprintf("max_restarts: %d\n", max_restarts);
fprintf("ada_max_restarts: %d\n", ada_max_restarts);
fprintf("truncation_length: %s\n", format_numeric_vector(truncation_length));
fprintf("max_num_quad_points: %d\n", max_num_quad_points);

fprintf("sketching_mat_type: %s\n", sketching_mat_type);
fprintf("sketching_size: %d\n", sketching_size);
fprintf("ada_sketcing_size_control: %d\n", ada_sketching_size_control);
fprintf("cond_tol: %.1e\n", cond_tol);

%% Load matrix.
load('./data/frac_laplacian/gnutella_comp.mat');
A = L;
N = size(A, 1);
b = A * b;
f_ex = ex_gnutella;

%% Choose parameters for the FUNM_QUAD restart algorithm.
addpath('funm_quad')
basic_param.function = 'invSqrt';
basic_param.restart_length = m;
basic_param.max_restarts = max_restarts;
basic_param.tol = quad_tol;
basic_param.transformation_parameter = 1;
basic_param.hermitian = 0;
basic_param.V_full = 0;
basic_param.H_full = 0;
basic_param.exact = f_ex;
basic_param.stopping_accuracy = stop_tol;
basic_param.inner_product = @(a,b) b'*a;
basic_param.thick = [];
basic_param.min_decay = 0.95;
basic_param.waitbar = 0;
basic_param.reorth_number = 0;
basic_param.truncation_length = inf;
basic_param.verbose = 1;
basic_param.check = check_flag;

add_param = construct_ada_param(...
    truncation_length, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    ada_max_restarts, ada_sketching_size_control, cond_tol);

%% Without implicit deflation.
result1 = run_single_method(A, b, f_ex, method_name, basic_param, add_param(1));
f1 = result1{1}.f;
out1 = result1{1}.out;
t1 = result1{1}.time;
fprintf("without implicit deflation, time: %.1e\n", t1);
fprintf("rel_err: %.1e\n", norm(f1 - f_ex) / norm(f_ex));
fprintf("num_cycle: %d\n", size(out1.appr, 2));
if check_flag == 1
    check_stats1 = get_max_check_stats(result1{1}.check_result);
    print_check_summary("without implicit deflation", check_stats1);
end

%% With implicit deflation.
basic_param.thick = @thick_quad;
basic_param.number_thick = 5;

result2 = run_single_method(A, b, f_ex, method_name, basic_param, add_param(1));
f2 = result2{1}.f;
out2 = result2{1}.out;
t2 = result2{1}.time;
fprintf("with implicit deflation, time: %.1e\n", t2);
fprintf("rel_err: %.1e\n", norm(f2 - f_ex) / norm(f_ex));
fprintf("num_cycle: %d\n", size(out2.appr, 2));
if check_flag == 1
    check_stats2 = get_max_check_stats(result2{1}.check_result);
    print_check_summary("with implicit deflation", check_stats2);
end

%% Plot
if truncation_length ~= inf
    displayname = method_name + " (t = " + string(truncation_length) + ")";
else
    displayname = method_name;
end

figure();
err1 = vecnorm(f_ex - out1.appr) / norm(f_ex);
semilogy(err1, 'g--+');
hold on;
for k = 2 : length(err1)
    text(k + 0.1, 2 * err1(k), num2str(out1.num_quadpoints(k)), ...
        'Color', [0 1 0], 'FontSize', 16, 'Rotation', 45);
end
err2 = vecnorm(f_ex - out2.appr) / norm(f_ex);
semilogy(err2, 'm--+');
for k = 2 : length(err2)
    text(k + 0.1, 2 * err2(k), num2str(out2.num_quadpoints(k)), ...
        'Color', [1 0 1], 'FontSize', 16, 'Rotation', 45);
end

legend('without deflation','implicit deflation')
xlabel('cycle'); ylabel('absolute 2-norm error');
title(displayname);
hold off

max_iter = max(size(out1.appr, 2), size(out2.appr, 2));
xticks(1 : ceil(max_iter / 10) : max_iter);

file_name = "implicit_deflation_rel_err_" + strrep(method_name, " ", "_") + "_" + string(m);
if truncation_length ~= inf
    file_name = file_name + "_" + string(truncation_length);
end
saveas(gcf, "./figure/frac_laplacian_invsqrt_thick_restart/" + file_name + ".eps", "epsc");

prefix_char = char(method_name);
ada_flag = strcmp(prefix_char(1:3), 'ada');
if ada_flag == 1
    figure();
    plot(out1.dim, '--o', "DisplayName", "without deflation");
    hold on;
    plot(out2.dim, '--p', "DisplayName", "implicit deflation");
    legend;
    xticks(1 : ceil(max_iter / 10) : max_iter);
    yticks(unique([out1.dim, out2.dim]));
    xlabel('cycle');
    ylabel('subspace dim');
    title(displayname);
    file_name = "implicit_deflation_subspace_dim_" + strrep(method_name, " ", "_") + "_" + string(m);
    if truncation_length ~= inf
        file_name = file_name + "_" + string(truncation_length);
    end
    saveas(gcf, "./figure/frac_laplacian_invsqrt_thick_restart/" + file_name + ".eps", "epsc");
end

function text = format_numeric_vector(values)
text = "[" + strjoin(string(values), ", ") + "]";
end

function print_check_summary(label, stats)
fprintf("  %s check max (before): AD=%s, Orth=%s, Cond=%s\n", ...
    label, ...
    format_check_value(stats.before_rel_err_AD), ...
    format_check_value(stats.before_rel_orth_err), ...
    format_check_value(stats.before_cond_V));
fprintf("  %s check max (after):  AD=%s, Orth=%s, Cond=%s\n", ...
    label, ...
    format_check_value(stats.after_rel_err_AD), ...
    format_check_value(stats.after_rel_orth_err), ...
    format_check_value(stats.after_cond_V));
end

function text = format_check_value(value)
if isempty(value)
    text = "[]";
else
    text = sprintf("%.1e", value);
end
end

function stats = get_max_check_stats(check_result)
stats = init_check_stats();
if isempty(check_result)
    return;
end

for cycle_idx = 1 : numel(check_result)
    cycle_result = check_result{cycle_idx};
    if isempty(cycle_result)
        continue;
    end
    stats = update_block_stats(stats, 'before', cycle_result.before);
    stats = update_block_stats(stats, 'after', cycle_result.after);
end
end

function stats = init_check_stats()
stats = struct();
stats.before_rel_err_AD = [];
stats.before_rel_orth_err = [];
stats.before_cond_V = [];
stats.after_rel_err_AD = [];
stats.after_rel_orth_err = [];
stats.after_cond_V = [];
end

function stats = update_block_stats(stats, prefix, block_result)
field_rel_err_AD = strcat(prefix, '_rel_err_AD');
field_rel_orth_err = strcat(prefix, '_rel_orth_err');
field_cond_V = strcat(prefix, '_cond_V');

stats.(field_rel_err_AD) = update_max_value(stats.(field_rel_err_AD), block_result.rel_err_AD);
stats.(field_rel_orth_err) = update_max_value(stats.(field_rel_orth_err), block_result.rel_orth_err);
stats.(field_cond_V) = update_max_value(stats.(field_cond_V), block_result.cond_V);
end

function max_value = update_max_value(max_value, value)
if isempty(value) || isnan(value)
    return;
end

value = abs(value);
if isempty(max_value)
    max_value = value;
else
    max_value = max(max_value, value);
end
end