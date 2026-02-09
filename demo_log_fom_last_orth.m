clear;
close all;
rng(2026);
maxNumCompThreads(1);

m = 50;
truncation_length = 5;
sk_type = "prod";
sk_factor = 2;
m_max = 100;
cond_tol = 1e6;

%% Initialize 2D Laplacian + some non-Herm part (no practical background).
nu = 1;
N = 500;
I = speye(N);
e = ones(N,1);
A = (N+1)^2*gallery('poisson',N);
s = eigs(A,1,'SM');
A = A/s;
o = ones(N,1);
D1 = (N+1)/2*spdiags([-o,0*o,o],-1:1,N,N);
D1 = kron(I,D1) + kron(D1,I);
A = A + nu * D1;
b = ones(N^2, 1);
b = b/norm(b);

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
[f_t, out_t] = funm_quad_fom_last_orth_tarnoldi(A,b,t_param);
t_t = toc;

num_it = length(out_t.num_quadpoints);
rel_err = norm(f - f_t) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_t);
t_rel_err0 = norm(out_t.appr(:, 1) - f) / norm(f);
fprintf("initial err: %e\n", t_rel_err0);
fprintf("\n\n");

fprintf("fom-s\n");
s_param = param;
s_param.sketch_dim_type = sk_type;
s_param.sketch_dim_factor = sk_factor;
tic;
[f_s, out_s] = funm_quad_fom_last_orth_sarnoldi(A,b,s_param);
t_s = toc;

num_it = length(out_s.num_quadpoints);
rel_err = norm(f - f_s) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_s);
s_rel_err0 = norm(out_s.appr(:, 1) - f) / norm(f);
fprintf("initial err: %e\n", s_rel_err0);
fprintf("\n\n");

fprintf("afom-t\n");
at_param = param;
at_param.truncation_length = truncation_length;
at_param.restart_length = m_max;
at_param.cond_tol = cond_tol;
tic;
[f_at, out_at] = funm_quad_ada_fom_last_orth_tarnoldi(A,b,at_param);
t_a = toc;

num_it = length(out_at.num_quadpoints);
rel_err = norm(f - f_at) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_a);
t_rel_err0 = norm(out_at.appr(:, 1) - f) / norm(f);
fprintf("initial err: %e\n", t_rel_err0);
fprintf("\n\n");

%% plot convergence curve and number of quadrature points
if ~isempty(out.appr)
    max_iter = max([length(out.appr), length(out_t.appr), length(out_s.appr)]);
    
    close all;
    figure();
    semilogy(vecnorm(f - out.appr) / norm(f), '--+', "DisplayName", "benchmark");
    hold on;
    semilogy(vecnorm(f - out_t.appr) / norm(f), '--x', "DisplayName", "fom-t");
    semilogy(vecnorm(f - out_s.appr) / norm(f), '--*', "DisplayName", "fom-s");
    semilogy(vecnorm(f - out_at.appr) / norm(f), '--o', "DisplayName", "afom-t");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('rel error compared to benchmark');

    figure();
    semilogy(out.update, '--+', "DisplayName", "benchmark");
    hold on;
    semilogy(out_t.update, '--x', "DisplayName", "fom-t");
    semilogy(out_s.update, '--*', "DisplayName", "fom-s");
    semilogy(out_at.update, '--o', "DisplayName", "afom-t");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('update norm');

    figure();
    plot(out.num_quadpoints, '--+', "DisplayName", "benchmark");
    hold on;
    plot(out_t.num_quadpoints, '--x', "DisplayName", "fom-t");
    plot(out_s.num_quadpoints, '--*', "DisplayName", "fom-s");
    plot(out_at.num_quadpoints, '--*', "DisplayName", "afom-t");
    legend;
    xticks(1 : max_iter)
    xlabel('cycle');
    ylabel('num of quad points');

    figure();
    plot(out_at.dim, '--*', "DisplayName", "afom-t");
    legend;
    xticks(1 : max_iter)
    xlabel('cycle');
    ylabel('num of subspace dim');
end
