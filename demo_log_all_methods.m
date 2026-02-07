clear;
close all;
rng(2026);
maxNumCompThreads(1);

m = 70;
truncation_length = 5;
sk_type = "prod";
sk_factor = 2;
ada_tol = inf;  % inf means always sketching
standard = "nonorth_fom";

%% Initialize 2D Laplacian + some non-Herm part (no practical background).
nu = 1e-12;
N = 500;
e = ones(N,1);
A = (N+1)^2*gallery('poisson',N);
I = speye(N);
o = ones(N,1);
D1 = (N+1)/2*spdiags([-o,0*o,o],-1:1,N,N);
D1 = kron(I,D1) + kron(D1,I);
A = A + nu * D1;
b = ones(N^2, 1);
b = b/norm(b);

load exact_solutions.mat;
f_ex = zeros(N^2, 1);
if nu == 0 && N == 40
    f_ex = exact_poisson_log;
end

%% choose parameters for the FUNM_QUAD restart algorithm
% jingyu: tol and stopping_accruacy are modified
addpath('funm_quad')
param.function = 'log';
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
param.reorth_number = 0;              % #reorthogonalizations
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
t_rel_err0 = norm(out_fom_t.appr(:, 1) - out.appr(:, 1)) / norm(out.appr(:, 1));
fprintf("initial err: %e\n", t_rel_err0);
fprintf("\n\n");

fprintf("sfom-t\n");
t_param = param;
t_param.truncation_length = truncation_length;
t_param.sketch_dim_type = sk_type;
t_param.sketch_dim_factor = sk_factor;
t_param.ada_tol = ada_tol;
t_param.standard = standard;
tic;
[f_sfom_t, out_sfom_t] = funm_quad_sfom_last_sorth_tarnoldi(A,b,t_param);
t_sfom_t = toc;

num_it = length(out_sfom_t.num_quadpoints);
rel_err = norm(f - f_sfom_t) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_sfom_t);
t_rel_err0 = norm(out_sfom_t.appr(:, 1) - out.appr(:, 1)) / norm(out.appr(:, 1));
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
s_rel_err0 = norm(out_fom_s.appr(:, 1) - out.appr(:, 1)) / norm(out.appr(:, 1));
fprintf("initial err: %e\n", s_rel_err0);
fprintf("\n\n");

fprintf("sfom-s\n");
s_param = param;
s_param.sketch_dim_type = sk_type;
s_param.sketch_dim_factor = sk_factor;
s_param.ada_tol = ada_tol;
s_param.standard = standard;
tic;
[f_sfom_s, out_sfom_s] = funm_quad_sfom_last_sorth_sarnoldi(A,b,s_param);
t_sfom_s = toc;

num_it = length(out_sfom_s.num_quadpoints);
rel_err = norm(f - f_sfom_s) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_sfom_s);
s_rel_err0 = norm(out_sfom_s.appr(:, 1) - out.appr(:, 1)) / norm(out.appr(:, 1));
fprintf("initial err: %e\n", s_rel_err0);
fprintf("number of sketching steps: %d\n", sum(out_sfom_s.sketching));
fprintf("\n\n");

%% print table
fprintf("\n\n");
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

%% plot convergence curve and number of quadrature points
if ~isempty(out.appr)
    max_iter = max([length(out.appr), ...
        length(out_fom_t.appr), ...
        length(out_fom_s.appr), ...
        length(out_sfom_t.appr), ...
        length(out_sfom_s.appr)]);
    
    close all;
    figure();
    semilogy(vecnorm(f - out.appr) / norm(f), 'g--+', "DisplayName", "benchmark");
    hold on;
    semilogy(vecnorm(f - out_fom_t.appr) / norm(f), 'r--x', "DisplayName", "fom-t");
    semilogy(vecnorm(f - out_fom_s.appr) / norm(f), 'b--*', "DisplayName", "fom-s");
    semilogy(vecnorm(f - out_sfom_t.appr) / norm(f), 'm--d', "DisplayName", "sfom-t");
    semilogy(vecnorm(f - out_sfom_s.appr) / norm(f), 'c--s', "DisplayName", "sfom-s");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('rel error compared to benchmark');

    figure();
    semilogy(out.update, 'g--+', "DisplayName", "benchmark");
    hold on;
    semilogy(out_fom_t.update, 'r--x', "DisplayName", "fom-t");
    semilogy(out_fom_s.update, 'b--*', "DisplayName", "fom-s");
    semilogy(out_sfom_t.update, 'm--d', "DisplayName", "sfom-t");
    semilogy(out_sfom_s.update, 'c--s', "DisplayName", "sfom-s");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('update norm');

    figure();
    plot(out.num_quadpoints, 'g--+', "DisplayName", "benchmark");
    hold on;
    plot(out_fom_t.num_quadpoints, 'r--x', "DisplayName", "fom-t");
    plot(out_fom_s.num_quadpoints, 'b--*', "DisplayName", "fom-s");
    plot(out_sfom_t.num_quadpoints, 'm--d', "DisplayName", "sfom-t");
    plot(out_sfom_s.num_quadpoints, 'c--s', "DisplayName", "sfom-s");
    legend;
    xticks(1 : max_iter)
    xlabel('cycle');
    ylabel('num of quad points');
end
