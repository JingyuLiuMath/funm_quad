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
truncation_length = [2, 1, 0];
max_num_quad_points = 8192;

sketching_mat_type = "sparse sign";
sketching_size = 2 * m;
ada_sketching_size_control = 2;
cond_tol = 1e4;

check_flag = 1;
num_test = 2;

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

%% choose parameters for the FUNM_QUAD restart algorithm
% jingyu: tol and stopping_accruacy are modified
addpath('funm_quad')
basic_param.function = 'invSqrt';
basic_param.restart_length = m;          % each restart cycle consists of 70 Arnoldi iterations
basic_param.max_restarts = max_restarts;            % perform at most 15 restart cycles
basic_param.tol = quad_tol;                   % tolerance for quadrature rule
basic_param.transformation_parameter = 1;     % parameter for the integral transformation
basic_param.hermitian = 0;                % the matrix A is Hermitian
basic_param.V_full = 0;                   % set 1 if you need Krylov basis
basic_param.H_full = 0;                   % do not store all Hessenberg matrices
basic_param.exact = [];
basic_param.stopping_accuracy = stop_tol;     % stopping accuracy
basic_param.inner_product = @(a,b) b'*a;  % use standard Euclidean inner product
basic_param.thick = [];                   % no implicit deflation is performed
basic_param.min_decay = .95;              % we desire linear error reduction of rate < .95 
basic_param.waitbar = 0;                  % show waitbar 
basic_param.reorth_number = 0;              % #reorthogonalizations
basic_param.truncation_length = inf;      % truncation length for Arnoldi 
basic_param.verbose = 1;                  % print information about progress of algorithm
basic_param.check = check_flag;

add_param = construct_ada_param(...
    truncation_length, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    ada_max_restarts, ada_sketching_size_control, cond_tol);

%% test
method_list = ["benchmark", ...
    "fix FOM-s", ...
    "ada FOM-t last-orth", "ada FOM-t last-sorth", "ada FOM-t whitening", ...
    "ada FOM-s last-orth", "ada FOM-s last-sorth", "ada FOM-s whitening", ...
    ];

num_methods = length(method_list);

% Initialize results storage
results_table = {};
results_table_names = ["Method", "Avg Error", "Bad Cases", "Time (s)", "Std", ...
    "B rel_err_AD", "B rel_orth_err", "B cond_V", ...
    "A rel_err_AD", "A rel_orth_err", "A cond_V"];
check_records = cell(num_methods, num_test);

% Test all methods
for method_idx = 1 : num_methods
    method_name = method_list(method_idx);
    
    fprintf("Testing method: %s\n", method_name);
    
    method_results = run_single_method(A, b, f_ex, method_name, basic_param, add_param);
    num_variants = numel(method_results);

    cnt = zeros(num_variants, 1);
    ava_err = zeros(num_variants, 1);
    ava_time = zeros(num_variants, 1);
    err_list = zeros(num_test, num_variants);
    max_check_stats = repmat(init_check_stats(), num_variants, 1);

    for it = 1 : num_test
        if it == 1
            current_results = method_results;
        else
            current_results = run_single_method(A, b, f_ex, method_name, basic_param, add_param);
        end

        check_records{method_idx, it} = current_results;

        for variant_idx = 1 : num_variants
            result = current_results{variant_idx};
            rel_err = result.rel_err;
            ava_err(variant_idx) = ava_err(variant_idx) + rel_err;
            ava_time(variant_idx) = ava_time(variant_idx) + result.time;
            err_list(it, variant_idx) = rel_err;

            if check_flag == 1 && isfield(result, "check_result")
                run_max_check_stats = get_max_check_stats(result.check_result);
                max_check_stats(variant_idx) = merge_check_stats(max_check_stats(variant_idx), run_max_check_stats);
            end

            if rel_err > 1e-4
                cnt(variant_idx) = cnt(variant_idx) + 1;
            end
        end
    end

    for variant_idx = 1 : num_variants
        avg_error = ava_err(variant_idx) / num_test;
        avg_time = ava_time(variant_idx) / num_test;
        std_error = std(err_list(:, variant_idx));
        result = method_results{variant_idx};
        method_disp = format_method_label(result.method, result.param.truncation_length);

        if has_valid_truncation_length(result)
            fprintf("  t = %g\n", result.param.truncation_length);
        end

        % Store results
        row_idx = size(results_table, 1) + 1;
        results_table{row_idx, 1} = method_disp;
        results_table{row_idx, 2} = avg_error;
        results_table{row_idx, 3} = cnt(variant_idx);
        results_table{row_idx, 4} = avg_time;
        results_table{row_idx, 5} = std_error;
        results_table{row_idx, 6} = max_check_stats(variant_idx).before_rel_err_AD;
        results_table{row_idx, 7} = max_check_stats(variant_idx).before_rel_orth_err;
        results_table{row_idx, 8} = max_check_stats(variant_idx).before_cond_V;
        results_table{row_idx, 9} = max_check_stats(variant_idx).after_rel_err_AD;
        results_table{row_idx, 10} = max_check_stats(variant_idx).after_rel_orth_err;
        results_table{row_idx, 11} = max_check_stats(variant_idx).after_cond_V;

        fprintf("  Avg Error: %.1e, Bad Cases: %d, Avg Time: %.1f s, Std: %.1e\n", ...
            avg_error, cnt(variant_idx), avg_time, std_error);
        fprintf("  Check max (before): AD=%s, Orth=%s, Cond=%s\n", ...
            format_check_value(max_check_stats(variant_idx).before_rel_err_AD), ...
            format_check_value(max_check_stats(variant_idx).before_rel_orth_err), ...
            format_check_value(max_check_stats(variant_idx).before_cond_V));
        fprintf("  Check max (after):  AD=%s, Orth=%s, Cond=%s\n", ...
            format_check_value(max_check_stats(variant_idx).after_rel_err_AD), ...
            format_check_value(max_check_stats(variant_idx).after_rel_orth_err), ...
            format_check_value(max_check_stats(variant_idx).after_cond_V));
        fprintf("\n");
    end
end

%% Print results table
fprintf("\n");
fprintf("%s\n", repmat('=', 1, 146));
fprintf("%-24s %8s %5s %8s %8s | %8s %8s %8s | %8s %8s %8s\n", ...
    "Method", "Avg", "Bad", "AvgTime", "Std", "B-AD", "B-Orth", "B-Cond", "A-AD", "A-Orth", "A-Cond");
fprintf("%s\n", repmat('-', 1, 146));

for method_idx = 1 : size(results_table, 1)
    method_disp = results_table{method_idx, 1};
    fprintf("%-36s %8.0e %5d %8.1f %8.0e | %8s %8s %8s | %8s %8s %8s\n", ...
        method_disp, ...
        results_table{method_idx, 2}, ...
        results_table{method_idx, 3}, ...
        results_table{method_idx, 4}, ...
        results_table{method_idx, 5}, ...
        format_check_value(results_table{method_idx, 6}), ...
        format_check_value(results_table{method_idx, 7}), ...
        format_check_value(results_table{method_idx, 8}), ...
        format_check_value(results_table{method_idx, 9}), ...
        format_check_value(results_table{method_idx, 10}), ...
        format_check_value(results_table{method_idx, 11}));
end

fprintf("%s\n", repmat('=', 1, 146));

%% Print LaTeX table (copy directly into .tex)
fprintf('\n');
fprintf('%%%% LaTeX table (part 1: summary)\n');
fprintf('\\begin{tabular}{lrrrr}\n');
fprintf('\\hline\n');
fprintf('Method & Avg & Bad & AvgTime & Std \\\\ \n');
fprintf('\\hline\n');
for method_idx = 1 : size(results_table, 1)
    method_tex = latex_escape(char(results_table{method_idx, 1}));
    fprintf('%s & %.1e & %d & %.1f & %.1e \\ \n', ...
        method_tex, ...
        results_table{method_idx, 2}, ...
        results_table{method_idx, 3}, ...
        results_table{method_idx, 4}, ...
        results_table{method_idx, 5});
end
fprintf('\\hline\n');
fprintf('\\end{tabular}\n');
fprintf('\n');
fprintf('%%%% LaTeX table (part 2: check)\n');
fprintf('\\begin{tabular}{lrrrrrr}\n');
fprintf('\\hline\n');
fprintf('Method & AD (before) & AD (after) & Orth (before) & Orth (after) & Cond (before) & Cond (after) \\\\ \n');
fprintf('\\hline\n');
for method_idx = 1 : size(results_table, 1)
    method_tex = latex_escape(char(results_table{method_idx, 1}));
    fprintf('%s & %s & %s & %s & %s & %s & %s \\\\ \n', ...
        method_tex, ...
        format_check_value(results_table{method_idx, 6}), ...
        format_check_value(results_table{method_idx, 9}), ...
        format_check_value(results_table{method_idx, 7}), ...
        format_check_value(results_table{method_idx, 10}), ...
        format_check_value(results_table{method_idx, 8}), ...
        format_check_value(results_table{method_idx, 11}));
end
fprintf('\\hline\n');
fprintf('\\end{tabular}\n');

function text = format_check_value(value)
if isempty(value)
    text = "[]";
else
    text = sprintf("%.1e", value);
end
end

function text = latex_escape(text)
text = strrep(text, '\\', '\\textbackslash{}');
text = strrep(text, '_', '\\_');
text = strrep(text, '%', '\\%');
text = strrep(text, '&', '\\&');
text = strrep(text, '#', '\\#');
end

function method_label = format_method_label(method_name, truncation_length)
method_label = char(method_name);
if ~isnumeric(truncation_length) || ~isscalar(truncation_length)
    return;
end

name_text = char(method_name);
if ~contains(name_text, "FOM-t")
    return;
end

method_label = sprintf('%s (t = %g)', name_text, truncation_length);
end

function text = format_numeric_vector(values)
text = "[" + strjoin(string(values), ", ") + "]";
end

function tf = has_valid_truncation_length(result)
tf = isfield(result, "param") && isfield(result.param, "truncation_length") && ...
    isnumeric(result.param.truncation_length) && isscalar(result.param.truncation_length) && ...
    isfinite(result.param.truncation_length);
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

function stats = merge_check_stats(stats, run_stats)
fields = fieldnames(stats);
for i = 1 : numel(fields)
    f = fields{i};
    stats.(f) = update_max_value(stats.(f), run_stats.(f));
end
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

