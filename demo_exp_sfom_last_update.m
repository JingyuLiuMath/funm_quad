clear;
close all;
rng(2026);
maxNumCompThreads(1);

truncation_length = 5;
sk_type = "prod";
sk_factor = 1.2;
ada_tol = inf;  % inf means always sketching
standard = "nonorth_fom";

%% Build discretization matrix for 2D convection-diffusion problem 
nu = 1;
N = 500;
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

%% choose parameters for the FUNM_QUAD restart algorithm
% jingyu: tol and stopping_accruacy are modified
addpath('funm_quad')
param.function = 'exp';
param.restart_length = 30;          % each restart cycle consists of 70 Arnoldi iterations
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
param.reorth_number = 1;            % #reorthogonalizations
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

fprintf("sfom with last update using t-arnoldi\n");
t_param = param;
t_param.truncation_length = truncation_length;
t_param.sketch_dim_type = sk_type;
t_param.sketch_dim_factor = sk_factor;
t_param.ada_tol = ada_tol;
t_param.standard = standard;
tic;
[f_t, out_t] = funm_quad_sfom_last_update_tarnoldi(A,b,t_param);
t_t = toc;

num_it = length(out_t.num_quadpoints);
rel_err = norm(f - f_t) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_t);
t_rel_err0 = norm(out_t.appr(:, 1) - out.appr(:, 1)) / norm(out.appr(:, 1));
fprintf("initial err: %e\n", t_rel_err0);
fprintf("number of sketching steps: %d\n", sum(out_t.sketching));
fprintf("\n\n");

fprintf("sfom with last update using s-arnoldi\n");
s_param = param;
s_param.sketch_dim_type = sk_type;
s_param.sketch_dim_factor = sk_factor;
s_param.ada_tol = ada_tol;
s_param.standard = standard;
tic;
[f_s, out_s] = funm_quad_sfom_last_update_sarnoldi(A,b,s_param);
t_s = toc;

num_it = length(out_s.num_quadpoints);
rel_err = norm(f - f_s) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t_s);
s_rel_err0 = norm(out_s.appr(:, 1) - out.appr(:, 1)) / norm(out.appr(:, 1));
fprintf("initial err: %e\n", s_rel_err0);
fprintf("number of sketching steps: %d\n", sum(out_s.sketching));
fprintf("\n\n");

%% plot convergence curve and number of quadrature points
if ~isempty(out.appr)
    max_iter = max([length(out.appr), length(out_t.appr), length(out_s.appr)]);
    
    close all;
    figure();
    semilogy(vecnorm(f - out.appr) / norm(f), 'g--+', "DisplayName", "benchmark");
    hold on;
    semilogy(vecnorm(f - out_t.appr) / norm(f), 'r--x', "DisplayName", "t-Arnoldi");
    semilogy(vecnorm(f - out_s.appr) / norm(f), 'b--*', "DisplayName", "s-Arnoldi");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('rel error compared to benchmark');

    figure();
    semilogy(out.update, 'g--+', "DisplayName", "benchmark");
    hold on;
    semilogy(out_t.update, 'r--x', "DisplayName", "t-Arnoldi");
    semilogy(out_s.update, 'b--*', "DisplayName", "s-Arnoldi");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('update norm');

    figure();
    plot(out.num_quadpoints, 'g--+', "DisplayName", "benchmark");
    hold on;
    plot(out_t.num_quadpoints, 'r--x', "DisplayName", "t-Arnoldi");
    plot(out_s.num_quadpoints, 'b--*', "DisplayName", "s-Arnoldi");
    legend;
    xticks(1 : max_iter)
    xlabel('cycle');
    ylabel('num of quad points');
end
