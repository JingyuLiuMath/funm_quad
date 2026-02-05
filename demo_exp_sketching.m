% Use the quadrature-based restart algorithm FUNM_QUAD described in
%
%  A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and 
%  stable Arnoldi restarts for matrix functions based on quadrature,
%  SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%
% to compute the exponential function of the 500x500 finite difference 
% discretization of a two-dimensional convection-diffusion problem for 
% convection parameters 0 (symmetric), 100 and 200 (non-symmetric) to 
% demonstrate the different behavior concerning speed of convergence 
% and required number of quadrature nodes for symmetric and non-symmetric 
% problems.

clear;
close all;
rng(2026);
maxNumCompThreads(1);

compare_mode = 0;
sk_type = "prod";
sk_factor = 4;
ada_tol = 1e-4;

%% Build discretization matrix for 2D convection-diffusion problem 
nu = 1;
N = 500;
D2 = (N+1)^2*gallery('tridiag',N);
I = speye(N);
D2 = kron(I,D2) + kron(D2,I);
o = ones(N,1);
D1 = (N+1)/2*spdiags([-o,0*o,o],-1:1,N,N);
D1 = kron(I,D1) + kron(D1,I);
A = D2 + nu*D1;

% choose time step s = 2*1e-3
s = 2*1e-3;
A = -s*A;

% choose right-hand side as normalized vector of all ones
b = ones(N^2,1); b = b/norm(b);


%% choose parameters for the FUNM_QUAD restart algorithm
% jingyu: tol and stopping_accruacy are modified
addpath('funm_quad')
param.function = 'exp';
param.restart_length = 70;          % each restart cycle consists of 70 Arnoldi iterations
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

fprintf("truncated Arnoldi\n");
sk1_param = param;
sk1_param.truncation_length = inf;
sk1_param.sketch_dim_type = sk_type;
sk1_param.sketch_dim_factor = sk_factor;
tic;
[sk1_f, sk1_out] = sketched_funm_quad_1(A,b,sk1_param);
sk1_t = toc;

num_it = length(sk1_out.num_quadpoints);
rel_err = norm(f - sk1_f) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, sk1_t);
sk1_rel_err0 = norm(sk1_out.appr(:, 1) - out.appr(:, 1)) / norm(out.appr(:, 1));
fprintf("initial err: %e\n", sk1_rel_err0);
fprintf("\n\n");

fprintf("adaptive truncated Arnoldi\n");
ada_sk1_param = param;
ada_sk1_param.truncation_length = 15;
ada_sk1_param.sketch_dim_type = sk_type;
ada_sk1_param.sketch_dim_factor = sk_factor;
ada_sk1_param.ada_tol = ada_tol;
tic;
[ada_sk1_f, ada_sk1_out] = sketched_funm_quad_1_adaptive(A,b,ada_sk1_param);
ada_sk1_t = toc;

num_it = length(ada_sk1_out.num_quadpoints);
rel_err = norm(f - ada_sk1_f) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, ada_sk1_t);
ada_sk1_rel_err0 = norm(ada_sk1_out.appr(:, 1) - out.appr(:, 1)) / norm(out.appr(:, 1));
fprintf("initial err: %e\n", ada_sk1_rel_err0);
fprintf("number of sketching steps: %d\n", sum(ada_sk1_out.sketching));
fprintf("\n\n");

fprintf("sketched Arnoldi\n");
sk2_param = param;
sk2_param.sketch_dim_type = sk_type;
sk2_param.sketch_dim_factor = sk_factor;
tic;
[sk2_f, sk2_out] = sketched_funm_quad_2(A,b,sk2_param);
sk2_t = toc;

num_it = length(sk2_out.num_quadpoints);
rel_err = norm(f - sk2_f) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, sk2_t);
sk2_rel_err0 = norm(sk2_out.appr(:, 1) - out.appr(:, 1)) / norm(out.appr(:, 1));
fprintf("initial err: %e\n", sk2_rel_err0);
fprintf("\n\n");

fprintf("adaptive sketched Arnoldi\n");
ada_sk2_param = param;
ada_sk2_param.sketch_dim_type = sk_type;
ada_sk2_param.sketch_dim_factor = sk_factor;
ada_sk2_param.ada_tol = ada_tol;
tic;
[ada_sk2_f, ada_sk2_out] = sketched_funm_quad_2_adaptive(A,b,ada_sk2_param);
sk2_t = toc;

num_it = length(ada_sk2_out.num_quadpoints);
rel_err = norm(f - ada_sk2_f) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, sk2_t);
ada_sk2_rel_err0 = norm(ada_sk2_out.appr(:, 1) - out.appr(:, 1)) / norm(out.appr(:, 1));
fprintf("initial err: %e\n", ada_sk2_rel_err0);
fprintf("number of sketching steps: %d\n", sum(ada_sk2_out.sketching));
fprintf("\n\n");

%% plot convergence curve and number of quadrature points
if ~isempty(out.appr)
    close all;
    figure();
    semilogy(vecnorm(f - out.appr) / norm(f), 'g--+', "DisplayName", "benchmark");
    hold on;
    semilogy(vecnorm(f - sk1_out.appr) / norm(f), 'r--x', "DisplayName", "t-Arnoldi");
    semilogy(vecnorm(f - sk2_out.appr) / norm(f), 'b--*', "DisplayName", "s-Arnoldi");
    semilogy(vecnorm(f - ada_sk1_out.appr) / norm(f), 'm--s', "DisplayName", "ada t-Arnoldi");
    semilogy(vecnorm(f - ada_sk2_out.appr) / norm(f), 'c--d', "DisplayName", "ada s-Arnoldi");
    legend;
    % xticks(1 : max([length(out.appr), length(sk1_out.appr), length(sk2_out.appr), length(ada_sk1_out.appr), length(ada_sk2_out.appr)]));
    xlabel('cycle');
    ylabel('rel error compared to benchmark');
    title('\nu = %d', nu);

    figure();
    semilogy(out.update, 'g--+', "DisplayName", "benchmark");
    hold on;
    semilogy(sk1_out.update, 'r--x', "DisplayName", "t-Arnoldi");
    semilogy(sk2_out.update, 'b--*', "DisplayName", "s-Arnoldi");
    semilogy(ada_sk1_out.update, 'm--s', "DisplayName", "ada t-Arnoldi");
    semilogy(ada_sk2_out.update, 'c--d', "DisplayName", "ada s-Arnoldi");
    legend;
    % xticks(1 : max([length(out.appr), length(sk1_out.appr), length(sk2_out.appr), length(ada_sk1_out.appr), length(ada_sk2_out.appr)]));
    xlabel('cycle');
    ylabel('update norm');
    title('\nu = %d', nu);

    figure();
    plot(out.num_quadpoints, 'g--+', "DisplayName", "benchmark");
    hold on;
    plot(sk1_out.num_quadpoints, 'r--x', "DisplayName", "t-Arnoldi");
    plot(sk2_out.num_quadpoints, 'b--*', "DisplayName", "s-Arnoldi");
    plot(ada_sk1_out.num_quadpoints, 'm--s', "DisplayName", "ada t-Arnoldi");
    plot(ada_sk2_out.num_quadpoints, 'c--d', "DisplayName", "ada s-Arnoldi");
    legend;
    % xticks(1 : max([length(out.appr), length(sk1_out.appr), length(sk2_out.appr), length(ada_sk1_out.appr), length(ada_sk2_out.appr)]));
    xlabel('cycle');
    ylabel('num of quad points');
    title('\nu = %d', nu);
end