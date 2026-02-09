clear;
close all;
rng(2026);
maxNumCompThreads(1);

% m = 70;
m = 200;
truncation_length = 5;
sk_type = "prod";
sk_factor = 2;
ada_tol = inf;  % inf means always sketching
standard = "nonorth_fom";
nu = 100;
N = 500;
m_max = 500;
cond_tol = 1e6;

fprintf("m: %d\n", m);
fprintf("truncation_length: %d\n", truncation_length);
fprintf("nu: %d\n", nu);
fprintf("N: %d\n", N);
fprintf("m_max: %d\n", m_max);
fprintf("cond_tol: %d\n", cond_tol);

%% Build discretization matrix for 2D convection-diffusion problem
D2 = (N+1)^2*gallery('tridiag',N);
I = speye(N);
D2 = kron(I,D2) + kron(D2,I);
o = ones(N,1);
D1 = (N+1)/2*spdiags([-o,0*o,o],-1:1,N,N);
D1 = kron(I,D1) + kron(D1,I);
A = D2 + nu * D1;

% choose time step s = 2*1e-3
s = 2*1e-3;
A = -s*A;

% choose right-hand side as normalized vector of all ones
b = ones(N^2, 1);
b = b/norm(b);

load exact_solutions.mat;
f_ex = zeros(N^2, 1);
if N == 500
    switch nu
        case 0
            f_ex = exact_convdiff_1;
        case 100
            f_ex = exact_convdiff_2;
        case 200
            f_ex = exact_convdiff_3;
    end
end

%% choose parameters for the FUNM_QUAD restart algorithm
% jingyu: tol and stopping_accruacy are modified
addpath('funm_quad')
param.function = 'exp';
param.restart_length = m;          % each restart cycle consists of 70 Arnoldi iterations
param.max_restarts = 15;            % perform at most 15 restart cycles
param.tol = 1e-7;                   % tolerance for quadrature rule
param.hermitian = 0;                % the matrix A is Hermitian
param.V_full = 0;                   % set 1 if you need Krylov basis
param.H_full = 0;                   % do not store all Hessenberg matrices
param.exact = [];
param.stopping_accuracy = 1e-8;     % stopping accuracy
param.inner_product = @(a,b) b'*a;  % use standard Euclidean inner product
param.thick = [];                   % no implicit deflation is performed
param.min_decay = .95;              % we desire linear error reduction of rate < .95
param.waitbar = 0;                  % show waitbar
param.reorth_number = 0;            % #reorthogonalizations
param.truncation_length = inf;      % truncation length for Arnoldi
param.verbose = 1;                  % print information about progress of algorithm

%% compute exp(A)b by quadrature-based restart algorithm
close all;
fprintf("\n\n");
fprintf("benchmark\n");
tic
[f,out] = funm_quad(A,b,param);
t = toc;

num_it = length(out.num_quadpoints);
rel_err = norm(f - f) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t)
fprintf("\n\n");

fprintf("fom-t\n");
t_param = param;
t_param.truncation_length = truncation_length;
tic;
[f_fom_t, out_fom_t] = funm_quad_fom_last_orth_tarnoldi(A,b,t_param);
t_fom_t = toc;

num_it = length(out_fom_t.num_quadpoints);
rel_err = norm(f - f_fom_t) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_fom_t);
t_rel_err0 = norm(out_fom_t.appr(:, 1) - f) / norm(f);
fprintf("initial err: %e\n", t_rel_err0);
fprintf("\n\n");

fprintf("sfom-t\n");
st_param = param;
st_param.truncation_length = truncation_length;
st_param.sketch_dim_type = sk_type;
st_param.sketch_dim_factor = sk_factor;
st_param.ada_tol = ada_tol;
st_param.standard = standard;
tic;
[f_sfom_t, out_sfom_t] = funm_quad_sfom_last_sorth_tarnoldi(A,b,st_param);
t_sfom_t = toc;

num_it = length(out_sfom_t.num_quadpoints);
rel_err = norm(f - f_sfom_t) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_sfom_t);
t_rel_err0 = norm(out_sfom_t.appr(:, 1) - f) / norm(f);
fprintf("initial err: %e\n", t_rel_err0);
fprintf("number of sketching steps: %d\n", sum(out_sfom_t.sketching));
fprintf("\n\n");

fprintf("fom-s\n");
s_param = param;
s_param.sketch_dim_type = sk_type;
s_param.sketch_dim_factor = sk_factor;
tic;
[f_fom_s, out_fom_s] = funm_quad_fom_last_orth_sarnoldi(A,b,s_param);
t_fom_s = toc;

num_it = length(out_fom_s.num_quadpoints);
rel_err = norm(f - f_fom_s) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_fom_s);
s_rel_err0 = norm(out_fom_s.appr(:, 1) - f) / norm(f);
fprintf("initial err: %e\n", s_rel_err0);
fprintf("\n\n");

fprintf("sfom-s\n");
ss_param = param;
ss_param.sketch_dim_type = sk_type;
ss_param.sketch_dim_factor = sk_factor;
ss_param.ada_tol = ada_tol;
ss_param.standard = standard;
tic;
[f_sfom_s, out_sfom_s] = funm_quad_sfom_last_sorth_sarnoldi(A,b,ss_param);
t_sfom_s = toc;

num_it = length(out_sfom_s.num_quadpoints);
rel_err = norm(f - f_sfom_s) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_sfom_s);
s_rel_err0 = norm(out_sfom_s.appr(:, 1) - f) / norm(f);
fprintf("initial err: %e\n", s_rel_err0);
fprintf("number of sketching steps: %d\n", sum(out_sfom_s.sketching));
fprintf("\n\n");

fprintf("afom-t\n");
at_param = param;
at_param.truncation_length = truncation_length;
at_param.restart_length = m_max;
at_param.cond_tol = 1e8;
tic;
[f_afom_t, out_afom_t] = funm_quad_ada_fom_last_orth_tarnoldi(A,b,at_param);
t_afom_t = toc;

num_it = length(out_afom_t.num_quadpoints);
rel_err = norm(f - f_afom_t) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_afom_t);
at_rel_err0 = norm(out_afom_t.appr(:, 1) - f) / norm(f);
fprintf("initial err: %e\n", at_rel_err0);
fprintf("\n\n");

%% save data
file_name = "./data/exp/exp_" + string(N) + "_" + string(m) + ".mat";
save(file_name, ...
    "f", "out", "t", ...
    "f_fom_t", "out_fom_t", "t_fom_t", ...
    "f_sfom_t", "out_sfom_t", "t_sfom_t", ...
    "f_fom_s", "out_fom_s", "t_fom_s", ...
    "f_sfom_s", "out_sfom_s", "t_sfom_s", ...
    "f_afom_t", "out_afom_t", "t_afom_t");

%% print table
fprintf("\n\n");

if norm(f_ex) ~= 0
    err_ex = norm(f - f_ex) / norm(f_ex);
    fprintf("err_ex: %e\n", err_ex);
end

num_it = length(out.num_quadpoints);
rel_err = norm(f - f) / norm(f);
fprintf("benchmark & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t);

num_it = length(out_fom_t.num_quadpoints);
rel_err = norm(f - f_fom_t) / norm(f);
fprintf("FOM-t & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t_fom_t);

num_it = length(out_sfom_t.num_quadpoints);
rel_err = norm(f - f_sfom_t) / norm(f);
fprintf("sFOM-t & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t_sfom_t);

num_it = length(out_fom_s.num_quadpoints);
rel_err = norm(f - f_fom_s) / norm(f);
fprintf("FOM-s & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t_fom_s);

num_it = length(out_sfom_s.num_quadpoints);
rel_err = norm(f - f_sfom_s) / norm(f);
fprintf("sFOM-s & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t_sfom_s);

num_it = length(out_afom_t.num_quadpoints);
rel_err = norm(f - f_afom_t) / norm(f);
fprintf("aFOM-t & %d & %.4e & %.4e \\\\ \n", num_it, rel_err, t_afom_t);

%% plot convergence curve and number of quadrature points
if ~isempty(out.appr)
    max_iter = max([length(out.appr), ...
        length(out_fom_t.appr), ...
        length(out_fom_s.appr), ...
        length(out_sfom_t.appr), ...
        length(out_sfom_s.appr), ...
        length(out_afom_t.appr)]);

    close all;
    figure();
    semilogy(vecnorm(f - out.appr) / norm(f), '--+', "DisplayName", "benchmark");
    hold on;
    semilogy(vecnorm(f - out_fom_t.appr) / norm(f), '--x', "DisplayName", "fom-t");
    semilogy(vecnorm(f - out_fom_s.appr) / norm(f), '--*', "DisplayName", "fom-s");
    semilogy(vecnorm(f - out_sfom_t.appr) / norm(f), '--d', "DisplayName", "sfom-t");
    semilogy(vecnorm(f - out_sfom_s.appr) / norm(f), '--s', "DisplayName", "sfom-s");
    semilogy(vecnorm(f - out_afom_t.appr) / norm(f), '--o', "DisplayName", "afom-t");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('rel error compared to benchmark');
    file_name = "exp_rel_err_" + string(N) + "_" + string(m);
    saveas(gcf, "./figure/exp/" + file_name + ".png", "png");
    saveas(gcf, "./figure/exp/" + file_name + ".eps", "epsc");


    figure();
    semilogy(out.update, '--+', "DisplayName", "benchmark");
    hold on;
    semilogy(out_fom_t.update, '--x', "DisplayName", "fom-t");
    semilogy(out_fom_s.update, '--*', "DisplayName", "fom-s");
    semilogy(out_sfom_t.update, '--d', "DisplayName", "sfom-t");
    semilogy(out_sfom_s.update, '--s', "DisplayName", "sfom-s");
    semilogy(out_afom_t.update, '--o', "DisplayName", "afom-t");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('update norm');
    file_name = "exp_norm_update_" + string(N) + "_" + string(m);
    saveas(gcf, "./figure/exp/" + file_name + ".png", "png");
    saveas(gcf, "./figure/exp/" + file_name + ".eps", "epsc");

    figure();
    plot(out.num_quadpoints, '--+', "DisplayName", "benchmark");
    hold on;
    plot(out_fom_t.num_quadpoints, '--x', "DisplayName", "fom-t");
    plot(out_fom_s.num_quadpoints, '--*', "DisplayName", "fom-s");
    plot(out_sfom_t.num_quadpoints, '--d', "DisplayName", "sfom-t");
    plot(out_sfom_s.num_quadpoints, '--s', "DisplayName", "sfom-s");
    plot(out_afom_t.num_quadpoints, '--o', "DisplayName", "afom-t");
    legend;
    xticks(1 : max_iter)
    xlabel('cycle');
    ylabel('num of quad points');
    file_name = "exp_num_quad_" + string(N) + "_" + string(m);
    saveas(gcf, "./figure/exp/" + file_name + ".png", "png");
    saveas(gcf, "./figure/exp/" + file_name + ".eps", "epsc");

    figure();
    plot(out_afom_t.dim, '--o', "DisplayName", "afom-t");
    xticks(1 : max_iter)
    xlabel('cycle');
    ylabel('num of subspace dim');
    file_name = "exp_subspace_dim_" + string(N) + "_" + string(m);
    saveas(gcf, "./figure/exp/" + file_name + ".png", "png");
    saveas(gcf, "./figure/exp/" + file_name + ".eps", "epsc");
end
