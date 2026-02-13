clear;
close all;
rng(2026);
maxNumCompThreads(1);

m = 300;
m_max = 600;
ada_max_restarts = 15;
truncation_length = 2;
max_num_quad_points = 1024;

sketching_mat_type = "Clarkson-Woodruff";
sketching_size = 2 * m;
ada_sketching_size_control = 2;
cond_tol = 1e6;

fprintf("m: %d\n", m);
fprintf("m_max: %d\n", m_max);
fprintf("ada_max_restarts: %d\n", ada_max_restarts);
fprintf("truncation_length: %d\n", truncation_length);
fprintf("max_num_quad_points: %d\n", max_num_quad_points);

fprintf("sketching_mat_type: %s\n", sketching_mat_type);
fprintf("sketching_size: %d\n", sketching_size);
fprintf("ada_sketcing_size_control: %d\n", ada_sketching_size_control);
fprintf("cond_tol: %d\n", cond_tol);

%% Load matrix.
load("./data/qcd/qcd_matrix_nonhermitian.mat");
N = size(Q, 1);
A = @(v) Q * (Q * v);
load("./data/qcd/qcd_nonhermitian_exact.mat");
c = zeros(N,1); c(1) = 1;
b = Q*c; normv = norm(b); b = b/norm(b);
f_ex = qcd_nonhermitian_exact / normv;

%% choose parameters for the FUNM_QUAD restart algorithm
% jingyu: tol and stopping_accruacy are modified
addpath('funm_quad')
param.function = 'invSqrt';
param.restart_length = m;          % each restart cycle consists of 70 Arnoldi iterations
param.max_restarts = 15;            % perform at most 15 restart cycles
param.tol = 1e-12;                   % tolerance for quadrature rule
param.transformation_parameter = 1;     % parameter for the integral transformation
param.hermitian = 0;                % the matrix A is Hermitian
param.V_full = 0;                   % set 1 if you need Krylov basis
param.H_full = 0;                   % do not store all Hessenberg matrices
param.exact = [];
param.stopping_accuracy = 1e-14;     % stopping accuracy
param.inner_product = @(a,b) b'*a;  % use standard Euclidean inner product
param.thick = [];                   % no implicit deflation is performed
param.min_decay = .95;              % we desire linear error reduction of rate < .95 
param.waitbar = 0;                  % show waitbar 
param.reorth_number = 0;              % #reorthogonalizations
param.truncation_length = inf;      % truncation length for Arnoldi 
param.verbose = 1;                  % print information about progress of algorithm

%% benchmark
close all;
fprintf("\n\n");
fprintf("benchmark\n");
tic
[f,out] = funm_quad(A,b,param);
t = toc;

num_it = size(out.appr, 2);
rel_err = norm(f_ex - f) / norm(f_ex);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t)
fprintf("\n\n");

%% fom-t
fprintf("fom-t\n");
fomt_param = param;
fomt_param.truncation_length = truncation_length;
fomt_param.max_num_quad_points = max_num_quad_points;
fomt_param.sarnoldi = 0;
fomt_param.last_update = "orth";
tic;
[f_fom_t, out_fom_t] = funm_quad_fom(A,b,fomt_param);
t_fom_t = toc;

num_it = size(out_fom_t.appr, 2);
rel_err = norm(f_ex - f_fom_t) / norm(f_ex);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_fom_t);
t_rel_err0 = norm(out_fom_t.appr(:, 1) - f_ex) / norm(f_ex);
fprintf("initial err: %e\n", t_rel_err0);
fprintf("\n\n");

%% sfom-t
fprintf("sfom-t\n");
sfomt_param = param;
sfomt_param.truncation_length = truncation_length;
sfomt_param.max_num_quad_points = max_num_quad_points;
sfomt_param.sketching_mat_type = sketching_mat_type;
sfomt_param.sketching_size = sketching_size;
sfomt_param.sarnoldi = 0;
sfomt_param.last_update = "sorth";
tic;
[f_sfom_t, out_sfom_t] = funm_quad_fom(A,b,sfomt_param);
t_sfom_t = toc;

num_it = size(out_sfom_t.appr, 2);
rel_err = norm(f_ex - f_sfom_t) / norm(f_ex);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_sfom_t);
t_rel_err0 = norm(out_sfom_t.appr(:, 1) - f_ex) / norm(f_ex);
fprintf("initial err: %e\n", t_rel_err0);
fprintf("\n\n");

%% fom-s
fprintf("fom-s\n");
foms_param = param;
foms_param.max_num_quad_points = max_num_quad_points;
foms_param.sketching_mat_type = sketching_mat_type;
foms_param.sketching_size = sketching_size;
foms_param.sarnoldi = 1;
foms_param.last_update = "orth";
tic;
[f_fom_s, out_fom_s] = funm_quad_fom(A,b,foms_param);
t_fom_s = toc;

num_it = size(out_fom_s.appr, 2);
rel_err = norm(f_ex - f_fom_s) / norm(f_ex);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_fom_s);
s_rel_err0 = norm(out_fom_s.appr(:, 1) - f_ex) / norm(f_ex);
fprintf("initial err: %e\n", s_rel_err0);
fprintf("\n\n");

%% sfom-s
fprintf("sfom-s\n");
sfoms_param = param;
sfoms_param.max_num_quad_points = max_num_quad_points;
sfoms_param.sketching_mat_type = sketching_mat_type;
sfoms_param.sketching_size = sketching_size;
sfoms_param.sarnoldi = 1;
sfoms_param.last_update = "sorth";
tic;
[f_sfom_s, out_sfom_s] = funm_quad_fom(A,b,sfoms_param);
t_sfom_s = toc;

num_it = size(out_sfom_s.appr, 2);
rel_err = norm(f_ex - f_sfom_s) / norm(f_ex);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_sfom_s);
s_rel_err0 = norm(out_sfom_s.appr(:, 1) - f_ex) / norm(f_ex);
fprintf("initial err: %e\n", s_rel_err0);
fprintf("\n\n");

%% afom-t
fprintf("afom-t\n");
afomt_param = param;
afomt_param.restart_length = m_max;
afomt_param.max_restarts = ada_max_restarts;
afomt_param.truncation_length = truncation_length;
afomt_param.max_num_quad_points = max_num_quad_points;
afomt_param.sketching_mat_type = sketching_mat_type;
afomt_param.ada_sketching_size_control = ada_sketching_size_control;
afomt_param.cond_tol = cond_tol;
afomt_param.last_update = "orth";
tic;
[f_afom_t, out_afom_t] = funm_quad_adaptive(A,b,afomt_param);
t_afom_t = toc;

num_it = size(out_afom_t.appr, 2);
rel_err = norm(f_ex - f_afom_t) / norm(f_ex);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_afom_t);
at_rel_err0 = norm(out_afom_t.appr(:, 1) - f_ex) / norm(f_ex);
fprintf("initial err: %e\n", at_rel_err0);
fprintf("\n\n");

%% asfom-t
fprintf("asfom-t\n");
asfomt_param = param;
asfomt_param.restart_length = m_max;
asfomt_param.max_restarts = ada_max_restarts;
asfomt_param.truncation_length = truncation_length;
asfomt_param.max_num_quad_points = max_num_quad_points;
asfomt_param.sketching_mat_type = sketching_mat_type;
asfomt_param.ada_sketching_size_control = ada_sketching_size_control;
asfomt_param.cond_tol = cond_tol;
asfomt_param.last_update = "sorth";
tic;
[f_asfom_t, out_asfom_t] = funm_quad_adaptive(A,b,asfomt_param);
t_asfom_t = toc;

num_it = size(out_asfom_t.appr, 2);
rel_err = norm(f_ex - f_asfom_t) / norm(f_ex);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_asfom_t);
at_rel_err0 = norm(out_asfom_t.appr(:, 1) - f_ex) / norm(f_ex);
fprintf("initial err: %e\n", at_rel_err0);
fprintf("\n\n");

%% save data
file_name = "./data/wikivote/wikivote_" + string(N) + "_" + string(m) + "_" + string(truncation_length) + ".mat";
save(file_name, ...
    "f", "out", "t", ...
    "f_fom_t", "out_fom_t", "t_fom_t", ...
    "f_sfom_t", "out_sfom_t", "t_sfom_t", ...
    "f_fom_s", "out_fom_s", "t_fom_s", ...
    "f_sfom_s", "out_sfom_s", "t_sfom_s", ...
    "f_afom_t", "out_afom_t", "t_afom_t", ...
    "f_asfom_t", "out_asfom_t", "t_asfom_t");

%% print table
fprintf("\n\n");

if norm(f_ex) ~= 0
    err_ex = norm(f - f_ex) / norm(f_ex);
    fprintf("err_ex: %e\n", err_ex);
end

num_it = size(out.appr, 2);
rel_err = norm(f_ex - f) / norm(f_ex);
fprintf("benchmark & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t);

num_it = size(out_fom_t.appr, 2);
rel_err = norm(f_ex - f_fom_t) / norm(f_ex);
fprintf("FOM-t & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t_fom_t);

num_it = size(out_sfom_t.appr, 2);
rel_err = norm(f_ex - f_sfom_t) / norm(f_ex);
fprintf("sFOM-t & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t_sfom_t);

num_it = size(out_fom_s.appr, 2);
rel_err = norm(f_ex - f_fom_s) / norm(f_ex);
fprintf("FOM-s & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t_fom_s);

num_it = size(out_sfom_s.appr, 2);
rel_err = norm(f_ex - f_sfom_s) / norm(f_ex);
fprintf("sFOM-s & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t_sfom_s);

num_it = size(out_afom_t.appr, 2);
rel_err = norm(f_ex - f_afom_t) / norm(f_ex);
fprintf("aFOM-t & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t_afom_t);

num_it = size(out_asfom_t.appr, 2);
rel_err = norm(f_ex - f_asfom_t) / norm(f_ex);
fprintf("asFOM-t & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t_asfom_t);

%% plot convergence curve and number of quadrature points
if ~isempty(out.appr)
    max_iter = max([size(out.appr, 2), ...
        size(out_fom_t.appr, 2), ...
        size(out_fom_s.appr, 2), ...
        size(out_sfom_t.appr, 2), ...
        size(out_sfom_s.appr, 2), ...
        size(out_afom_t.appr, 2), ...
        size(out_asfom_t.appr, 2)]);

    close all;
    figure();
    semilogy(vecnorm(f_ex - out.appr) / norm(f_ex), '--+', "DisplayName", "benchmark");
    hold on;
    semilogy(vecnorm(f_ex - out_fom_t.appr) / norm(f_ex), '--x', "DisplayName", "fom-t");
    semilogy(vecnorm(f_ex - out_fom_s.appr) / norm(f_ex), '--*', "DisplayName", "fom-s");
    semilogy(vecnorm(f_ex - out_sfom_t.appr) / norm(f_ex), '--d', "DisplayName", "sfom-t");
    semilogy(vecnorm(f_ex - out_sfom_s.appr) / norm(f_ex), '--s', "DisplayName", "sfom-s");
    semilogy(vecnorm(f_ex - out_afom_t.appr) / norm(f_ex), '--o', "DisplayName", "afom-t");
    semilogy(vecnorm(f_ex - out_asfom_t.appr) / norm(f_ex), '--p', "DisplayName", "asfom-t");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('rel error compared to benchmark');
    file_name = "qcd_rel_err_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/qcd/" + file_name + ".png", "png");
    saveas(gcf, "./figure/qcd/" + file_name + ".eps", "epsc");


    figure();
    semilogy(out.update, '--+', "DisplayName", "benchmark");
    hold on;
    semilogy(out_fom_t.update, '--x', "DisplayName", "fom-t");
    semilogy(out_fom_s.update, '--*', "DisplayName", "fom-s");
    semilogy(out_sfom_t.update, '--d', "DisplayName", "sfom-t");
    semilogy(out_sfom_s.update, '--s', "DisplayName", "sfom-s");
    semilogy(out_afom_t.update, '--o', "DisplayName", "afom-t");
    semilogy(out_asfom_t.update, '--p', "DisplayName", "asfom-t");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('update norm');
    file_name = "qcd_norm_update_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/qcd/" + file_name + ".png", "png");
    saveas(gcf, "./figure/qcd/" + file_name + ".eps", "epsc");

    figure();
    plot(out.num_quadpoints, '--+', "DisplayName", "benchmark");
    hold on;
    plot(out_fom_t.num_quadpoints, '--x', "DisplayName", "fom-t");
    plot(out_fom_s.num_quadpoints, '--*', "DisplayName", "fom-s");
    plot(out_sfom_t.num_quadpoints, '--d', "DisplayName", "sfom-t");
    plot(out_sfom_s.num_quadpoints, '--s', "DisplayName", "sfom-s");
    plot(out_afom_t.num_quadpoints, '--o', "DisplayName", "afom-t");
    plot(out_asfom_t.num_quadpoints, '--p', "DisplayName", "asfom-t");
    legend;
    xticks(1 : max_iter)
    xlabel('cycle');
    ylabel('num of quad points');
    file_name = "qcd_num_quad_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/qcd/" + file_name + ".png", "png");
    saveas(gcf, "./figure/qcd/" + file_name + ".eps", "epsc");

    figure();
    plot(out_afom_t.dim, '--o', "DisplayName", "afom-t");
    hold on;
    plot(out_asfom_t.dim, '--p', "DisplayName", "asfom-t");
    legend;
    xticks(1 : max_iter)
    xlabel('cycle');
    ylabel('num of subspace dim');
    file_name = "qcd_subspace_dim_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/qcd/" + file_name + ".png", "png");
    saveas(gcf, "./figure/qcd/" + file_name + ".eps", "epsc");
end