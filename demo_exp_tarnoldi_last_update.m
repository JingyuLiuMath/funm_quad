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
sk_factor = 1.2;
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
A = D2 + nu * D1;

% choose time step s = 2*1e-3
s = 2*1e-3;
A = -s*A;

% choose right-hand side as normalized vector of all ones
b = ones(N^2,1); b = b/norm(b);


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
if compare_mode == 1
    param.exact = exact_convdiff_1;     % exact solution. If not known set to []
else
    param.exact = [];
end
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

fprintf("truncated Arnoldi with last update\n");
new_param = param;
new_param.truncation_length = 10;
tic;
[f1, out1] = funm_quad_last_update(A,b,new_param);
t1 = toc;

num_it = length(out1.num_quadpoints);
rel_err = norm(f - f1) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t1);
sk1_rel_err0 = norm(out1.appr(:, 1) - out.appr(:, 1)) / norm(out.appr(:, 1));
fprintf("initial err: %e\n", sk1_rel_err0);
fprintf("\n\n");


%% plot convergence curve and number of quadrature points
max_iter = max([length(out.appr), length(out1.appr)]);
if ~isempty(out.appr)
    close all;
    figure();
    semilogy(vecnorm(f - out.appr) / norm(f), 'g--+', "DisplayName", "benchmark");
    hold on;
    semilogy(vecnorm(f - out1.appr) / norm(f), 'r--x', "DisplayName", "last updated t-Arnoldi");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('rel error compared to benchmark');

    figure();
    semilogy(out.update, 'g--+', "DisplayName", "benchmark");
    hold on;
    semilogy(out1.update, 'r--x', "DisplayName", "last updated t-Arnoldi");
    legend;
    xticks(1 : max_iter);
    xlabel('cycle');
    ylabel('update norm');

    figure();
    plot(out.num_quadpoints, 'g--+', "DisplayName", "benchmark");
    hold on;
    plot(out1.num_quadpoints, 'r--x', "DisplayName", "last updated t-Arnoldi");
    legend;
    xticks(1 : max_iter)
    xlabel('cycle');
    ylabel('num of quad points');
end
