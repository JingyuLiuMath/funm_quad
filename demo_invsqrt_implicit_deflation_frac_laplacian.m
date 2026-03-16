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

method_list = ["benchmark", ...
    "fix FOM-t last-orth", "fix FOM-t last-sorth", ...
    "fix FOM-s", ...
    "ada FOM-t last-orth", "ada FOM-t last-sorth"];
method_name = "fix FOM-t last-orth";

%% Load matrix.
load('./data/frac_laplacian/gnutella_comp.mat');
A = L;
N = size(A, 1);
b = A * b;
f_ex = ex_gnutella;

%% Choose parameters for the FUNM_QUAD algorithm (no implicit deflation)
addpath('funm_quad')
basic_param.function = 'invSqrt';
basic_param.restart_length = m;              % Each restart cycle consists of 20 Arnoldi iterations
basic_param.max_restarts = max_restarts;                % We perform at most 20 restart cycles
basic_param.tol = quad_tol;                      % Tolerance for quadrature rule
basic_param.transformation_parameter = 1;
basic_param.hermitian = 0;                    % set 0 if A is not Hermitian
basic_param.V_full = 0;                       % set 1 if you need Krylov basis
basic_param.H_full = 0;                       % do not store all Hessenberg matrices
basic_param.exact = [];        % Exact solution. If not known set to []
basic_param.stopping_accuracy = stop_tol;        % stopping accuracy
basic_param.inner_product = @(a,b) b'*a;      % Use standard euclidean inner product
basic_param.thick = [];                       % No implicit deflation is performed
basic_param.min_decay = 0.95;                 % we desire linear error reduction of rate < .95
basic_param.waitbar = 0;                      % show waitbar
basic_param.reorth_number = 0;                % #reorthogonalizations
basic_param.truncation_length = inf;          % truncation length for Arnoldi
basic_param.verbose = 1;                      % print information about progress of algorithm

add_param = construct_ada_param(...
    truncation_length, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    ada_max_restarts, ada_sketching_size_control, cond_tol);

%% without implicit deflation
result1 = run_single_method(A, b, f_ex, method_name, basic_param, add_param);

f1 = result1{1}.f;
out1 = result1{1}.out;
t1 = result1{1}.time;

%% with implicit deflation
basic_param.thick = @thick_quad;              % Thick restart function for implicit deflation
basic_param.number_thick = 5;                 % Number of target eigenvalues for implicit deflation
% basic_param.restart_length = m - 5;

result2 = run_single_method(A, b, f_ex, method_name, basic_param, add_param);
f2 = result2{1}.f;
out2 = result2{1}.out;
t2 = result2{1}.time;

%% Print
num_it = size(out1.appr, 2);
rel_err = norm(f_ex - f1) / norm(f_ex);
fprintf("%s (without deflation) & %d & %.4e & %.4e \\\\ \n", method_name, num_it, rel_err, t1);

num_it = size(out2.appr, 2);
rel_err = norm(f_ex - f2) / norm(f_ex);
fprintf("%s (implicit deflation) & %d & %.4e & %.4e \\\\ \n", method_name, num_it, rel_err, t2);

%% plot
max_iter = max(size(out1.appr, 2), size(out2.appr, 2));

figure();
semilogy(vecnorm(f_ex - out1.appr) / norm(f_ex), '--o', "DisplayName", "without deflation");
hold on;
semilogy(vecnorm(f_ex - out2.appr) / norm(f_ex), '--p', "DisplayName", "implicit deflation");
legend;
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('rel error compared to exact');
% if method == "benchmark"
%     file_name = "frac_laplacian_id_rel_err_" + method + "_" + string(N) + "_" + string(m);
%     saveas(gcf, "./figure/frac_laplacian_implicit_deflation/" + file_name + ".eps", "epsc");
% else
%     file_name = "frac_laplacian_id_rel_err_" + method + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
%     saveas(gcf, "./figure/frac_laplacian_implicit_deflation/" + file_name + ".eps", "epsc");
% 
% end

figure();
semilogy(out1.update, '--o', "DisplayName", "without deflation");
hold on;
semilogy(out2.update, '--p', "DisplayName", "implicit deflation");
legend;
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('update norm');
% if method == "benchmark"
%     file_name = "frac_laplacian_id_norm_update_" + method + "_" + string(N) + "_" + string(m);
%     saveas(gcf, "./figure/frac_laplacian_implicit_deflation/" + file_name + ".eps", "epsc");
% else
%     file_name = "frac_laplacian_id_norm_update_" + method + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
%     saveas(gcf, "./figure/frac_laplacian_implicit_deflation/" + file_name + ".eps", "epsc");
% 
% end

figure();
plot(out1.num_quadpoints, '--o', "DisplayName", "without deflation");
hold on;
plot(out2.num_quadpoints, '--p', "DisplayName", "implicit deflation");
legend;
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('num of quad points');
% if method == "benchmark"
%     file_name = "frac_laplacian_id_num_quad_" + method + "_" + string(N) + "_" + string(m);
%     saveas(gcf, "./figure/frac_laplacian_implicit_deflation/" + file_name + ".eps", "epsc");
% else
%     file_name = "frac_laplacian_id_num_quad_" + method + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
%     saveas(gcf, "./figure/frac_laplacian_implicit_deflation/" + file_name + ".eps", "epsc");
% 
% end

prefix_char = char(method_name);
ada_flag = strcmp(prefix_char(1:3), 'ada');

if ada_flag
    figure();
    plot(out1.dim, '--o', "DisplayName", "without deflation");
    hold on;
    plot(out2.dim, '--p', "DisplayName", "implicit deflation");
    legend;
    xticks(1 : ceil(max_iter / 10) : max_iter);
    yticks(unique([out1.dim, out2.dim]));
    xlabel('cycle');
    ylabel('subspace dim');
    file_name = "frac_laplacian_id_subspace_dim_" + method + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/frac_laplacian_implicit_deflation/" + file_name + ".eps", "epsc");
end