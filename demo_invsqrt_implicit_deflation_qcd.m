clear;
close all;
rng(2026);
maxNumCompThreads(1);
warning off;

quad_tol = 1e-7;
stop_tol = 1e-8;

m = 30;
max_restarts = 15;
ada_max_restarts = 200;
truncation_length = 1;
max_num_quad_points = 8192;

sketching_mat_type = "sparse sign";
sketching_size = 2 * m;
ada_sketching_size_control = 2;
cond_tol = 1e6;

method_num = 3;
switch method_num
    case 1
        method = "benchmark";
    case 2
        method = "aFOM-t";
    case 3
        method = "asFOM-t";
end

%% Load matrix.
load("./data/qcd/qcd_matrix_nonhermitian.mat");
N = size(Q, 1);
A = @(v) Q * (Q * v);
load("./data/qcd/qcd_nonhermitian_exact.mat");
c = zeros(N,1); c(1) = 1;
b = Q*c; normv = norm(b); b = b/norm(b);
f_ex = qcd_nonhermitian_exact / normv;

%% Choose parameters for the FUNM_QUAD algorithm (no implicit deflation)
addpath('funm_quad')
param.function = 'invSqrt';
param.restart_length = m;              % Each restart cycle consists of 20 Arnoldi iterations
param.max_restarts = max_restarts;                % We perform at most 20 restart cycles
param.tol = quad_tol;                      % Tolerance for quadrature rule
param.transformation_parameter = 1;
param.hermitian = 0;                    % set 0 if A is not Hermitian
param.V_full = 0;                       % set 1 if you need Krylov basis
param.H_full = 0;                       % do not store all Hessenberg matrices
param.exact = [];        % Exact solution. If not known set to []
param.stopping_accuracy = stop_tol;        % stopping accuracy
param.inner_product = @(a,b) b'*a;      % Use standard euclidean inner product
param.thick = [];                       % No implicit deflation is performed
param.min_decay = 0.95;                 % we desire linear error reduction of rate < .95
param.waitbar = 0;                      % show waitbar
param.reorth_number = 0;                % #reorthogonalizations
param.truncation_length = inf;          % truncation length for Arnoldi
param.verbose = 1;                      % print information about progress of algorithm

%% without implicit deflation
if method == "benchmark"
    tic
    [f1,out1] = funm_quad(A,b,param);
    t1 = toc;
elseif method == "aFOM-t"
    afomt_param = param;
    afomt_param.restart_length = m;
    afomt_param.max_restarts = ada_max_restarts;
    afomt_param.truncation_length = truncation_length;
    afomt_param.max_num_quad_points = max_num_quad_points;
    afomt_param.sketching_mat_type = sketching_mat_type;
    afomt_param.ada_sketching_size_control = ada_sketching_size_control;
    afomt_param.cond_tol = cond_tol;
    afomt_param.last_update = "orth";

    tic
    [f1,out1] = funm_quad_adaptive(A,b,afomt_param);
    t1 = toc;
elseif method == "asFOM-t"
    asfomt_param = param;
    asfomt_param.restart_length = m;
    asfomt_param.max_restarts = ada_max_restarts;
    asfomt_param.truncation_length = truncation_length;
    asfomt_param.max_num_quad_points = max_num_quad_points;
    asfomt_param.sketching_mat_type = sketching_mat_type;
    asfomt_param.ada_sketching_size_control = ada_sketching_size_control;
    asfomt_param.cond_tol = cond_tol;
    asfomt_param.last_update = "sorth";

    tic
    [f1,out1] = funm_quad_adaptive(A,b,asfomt_param);
    t1 = toc;
end

%% with implicit deflation
param.thick = @thick_quad;              % Thick restart function for implicit deflation
param.number_thick = 5;                 % Number of target eigenvalues for implicit deflation

if method == "benchmark"
    tic
    [f2,out2] = funm_quad(A,b,param);
    t2 = toc;
elseif method == "aFOM-t"
    afomt_param = param;
    afomt_param.restart_length = m;
    afomt_param.max_restarts = ada_max_restarts;
    afomt_param.truncation_length = truncation_length;
    afomt_param.max_num_quad_points = max_num_quad_points;
    afomt_param.sketching_mat_type = sketching_mat_type;
    afomt_param.ada_sketching_size_control = ada_sketching_size_control;
    afomt_param.cond_tol = cond_tol;
    afomt_param.last_update = "orth";

    tic
    [f2,out2] = funm_quad_adaptive(A,b,afomt_param);
    t2 = toc;
elseif method == "asFOM-t"
    asfomt_param = param;
    asfomt_param.restart_length = m;
    asfomt_param.max_restarts = ada_max_restarts;
    asfomt_param.truncation_length = truncation_length;
    asfomt_param.max_num_quad_points = max_num_quad_points;
    asfomt_param.sketching_mat_type = sketching_mat_type;
    asfomt_param.ada_sketching_size_control = ada_sketching_size_control;
    asfomt_param.cond_tol = cond_tol;
    asfomt_param.last_update = "sorth";

    tic
    [f2,out2] = funm_quad_adaptive(A,b,asfomt_param);
    t2 = toc;
end

%% Print
num_it = size(out1.appr, 2);
rel_err = norm(f_ex - f1) / norm(f_ex);
fprintf("%s (without deflation) & %d & %.4e & %.4e \\\\ \n", method, num_it, rel_err, t1);

num_it = size(out2.appr, 2);
rel_err = norm(f_ex - f2) / norm(f_ex);
fprintf("%s (implicit deflation) & %d & %.4e & %.4e \\\\ \n", method, num_it, rel_err, t2);

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
if method == "benchmark"
    file_name = "qcd_id_rel_err_" + method + "_" + string(N) + "_" + string(m);
    saveas(gcf, "./figure/qcd_implicit_deflation/" + file_name + ".eps", "epsc");
else
    file_name = "qcd_id_rel_err_" + method + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/qcd_implicit_deflation/" + file_name + ".eps", "epsc");

end

figure();
semilogy(out1.update, '--o', "DisplayName", "without deflation");
hold on;
semilogy(out2.update, '--p', "DisplayName", "implicit deflation");
legend;
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('update norm');
if method == "benchmark"
    file_name = "qcd_id_norm_update_" + method + "_" + string(N) + "_" + string(m);
    saveas(gcf, "./figure/qcd_implicit_deflation/" + file_name + ".eps", "epsc");
else
    file_name = "qcd_id_norm_update_" + method + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/qcd_implicit_deflation/" + file_name + ".eps", "epsc");

end

figure();
plot(out1.num_quadpoints, '--o', "DisplayName", "without deflation");
hold on;
plot(out2.num_quadpoints, '--p', "DisplayName", "implicit deflation");
legend;
xticks(1 : ceil(max_iter / 10) : max_iter);
xlabel('cycle');
ylabel('num of quad points');
if method == "benchmark"
    file_name = "qcd_id_num_quad_" + method + "_" + string(N) + "_" + string(m);
    saveas(gcf, "./figure/qcd_implicit_deflation/" + file_name + ".eps", "epsc");
else
    file_name = "qcd_id_num_quad_" + method + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/qcd_implicit_deflation/" + file_name + ".eps", "epsc");

end

if method ~= "benchmark"
    figure();
    plot(out1.dim, '--o', "DisplayName", "without deflation");
    hold on;
    plot(out2.dim, '--p', "DisplayName", "implicit deflation");
    legend;
    xticks(1 : ceil(max_iter / 10) : max_iter);
    yticks(unique([out1.dim, out2.dim]));
    xlabel('cycle');
    ylabel('subspace dim');
    file_name = "qcd_id_subspace_dim_" + method + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/qcd_implicit_deflation/" + file_name + ".eps", "epsc");
end