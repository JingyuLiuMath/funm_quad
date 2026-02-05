clear;
close all;
rng(2026);
maxNumCompThreads(1);

sk_type = "prod";
sk_factor = 1.2;

%% Build discretization matrix for 2D convection-diffusion problem 
N = 500;
D2 = gallery('tridiag',N);
I = speye(N);
D2 = kron(I,D2) + kron(D2,I);
A = sprandn(D2);
herm_err = norm(A - A', "fro")/ norm(A, "fro");
if herm_err == 0
    fprintf("Herm mat!\n");
end
    
% choose right-hand side as normalized vector of all ones
b = ones(N^2,1); b = b/norm(b);

%% choose parameters for the FUNM_QUAD restart algorithm
addpath('funm_quad')
param.function = 'exp';
param.restart_length = 70;          % each restart cycle consists of 70 Arnoldi iterations
param.max_restarts = 15;            % perform at most 15 restart cycles
param.tol = 1e-7;                  % tolerance for quadrature rule
param.hermitian = 0;                % the matrix A is Hermitian
param.V_full = 0;                   % set 1 if you need Krylov basis
param.H_full = 0;                   % do not store all Hessenberg matrices
param.exact = [];     % exact solution. If not known set to []
param.stopping_accuracy = 1e-8;    % stopping accuracy
param.inner_product = @(a,b) b'*a;  % use standard Euclidean inner product
param.thick = [];                   % no implicit deflation is performed
param.min_decay = .95;              % we desire linear error reduction of rate < .95 
param.waitbar = 0;                  % show waitbar 
param.reorth_number = 0;            % #reorthogonalizations
param.truncation_length = inf;      % truncation length for Arnoldi 
param.verbose = 1;                  % print information about progress of algorithm

%% compute exp(A)b by quadrature-based restart algorithm
fprintf("\n\n");

fprintf("benchmark\n");
tic
[f, out] = funm_quad(A,b,param);
t = toc;

num_it = length(out.num_quadpoints);
rel_err = norm(f - f) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, t)
fprintf("\n\n");

fprintf("truncated Arnoldi\n");
sk1_param = param;
sk1_param.truncation_length = 70;
sk1_param.sketch_dim_type = sk_type;
sk1_param.sketch_dim_factor = sk_factor;
tic;
[sk1_f, sk1_out] = sketched_funm_quad_1(A,b,sk1_param);
sk1_t = toc;

num_it = length(sk1_out.num_quadpoints);
rel_err = norm(f - sk1_f) / norm(f);
fprintf("iter rel_err time\n");
fprintf(" %d & %.4e & %.4e \n", num_it, rel_err, sk1_t);
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
fprintf("\n\n");

%%
if ~isempty(out.appr)
    figure();
    semilogy(vecnorm(f - out.appr) / norm(f), 'g--+', "DisplayName", "benchmark");
    hold on;
    semilogy(vecnorm(f - sk1_out.appr) / norm(f), 'r--x', "DisplayName", "t-Arnoldi");
    semilogy(vecnorm(f - sk2_out.appr) / norm(f), 'b--*', "DisplayName", "s-GS");
    legend;
    xticks(1 : max([length(out.appr), length(sk1_out.appr), length(sk2_out.appr)]))
    xlabel('cycle');
    ylabel('rel error compared to benchmark');

    figure();
    semilogy(out.update, 'g--+', "DisplayName", "benchmark");
    hold on;
    semilogy(sk1_out.update, 'r--x', "DisplayName", "t-Arnoldi");
    semilogy(sk2_out.update, 'b--*', "DisplayName", "s-GS");
    legend;
    xticks(1 : max([length(out.appr), length(sk1_out.appr), length(sk2_out.appr)]))
    xlabel('cycle');
    ylabel('update norm');

    figure();
    plot(out.num_quadpoints, 'g--+', "DisplayName", "benchmark");
    hold on;
    plot(sk1_out.num_quadpoints, 'r--x', "DisplayName", "t-Arnoldi");
    plot(sk2_out.num_quadpoints, 'b--*', "DisplayName", "s-GS");
    legend;
    xticks(1 : max([length(out.appr), length(sk1_out.appr), length(sk2_out.appr)]))
    xlabel('cycle');
    ylabel('num of quad points');
end
