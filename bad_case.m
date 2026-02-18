clear;
close all;
rng(2026);
maxNumCompThreads(1);

% m = 50;
% m_max = 100;
m = 100;
m_max = 200;
ada_max_restarts = 15;
truncation_length = 0;
max_num_quad_points = 100000;

sketching_mat_type = "Clarkson-Woodruff";
sketching_size = 4 * m;
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
load('data/wikivote/wiki-Vote.mat');
load('data/wikivote/wiki-Vote-comp.mat');
A = -Problem.A; ee = -ee; N = size(A,1); 
b = ones(N,1);
f_ex = ex_expm;

%% choose parameters for the FUNM_QUAD restart algorithm
% jingyu: tol and stopping_accruacy are modified
addpath('funm_quad')
param.function = 'exp';
param.restart_length = m;          % each restart cycle consists of 70 Arnoldi iterations
param.max_restarts = 15;            % perform at most 15 restart cycles
param.tol = 1e-10;                   % tolerance for quadrature rule
param.hermitian = 0;                % the matrix A is Hermitian
param.V_full = 0;                   % set 1 if you need Krylov basis
param.H_full = 0;                   % do not store all Hessenberg matrices
param.exact = [];
param.stopping_accuracy = 1e-14;     % stopping accuracy
param.inner_product = @(a,b) b'*a;  % use standard Euclidean inner product
param.thick = [];                   % no implicit deflation is performed
param.min_decay = .95;              % we desire linear error reduction of rate < .95
param.waitbar = 0;                  % show waitbar
param.reorth_number = 0;            % #reorthogonalizations
param.truncation_length = inf;      % truncation length for Arnoldi
param.verbose = 2;                  % print information about progress of algorithm

%% afom-t
% rng(20260213)  % quadrature points were not enough in iter 4, bad acc in iter 4.
% rng(2026)  % quadrature points were not enough in iter 4, bad acc in iter 4.
% rng(1)  % quadrature points were not enough in iter 4, bad acc in iter 4.
% rng(1024)  % good case.
% rng(0524) % quadrature points were not enough in iter 4, bad acc in iter 4.

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
iter_err = vecnorm(f_ex - out_afom_t.appr) / norm(f_ex);
fprintf("rel err in each iter\n");
for j = 1 : length(iter_err)
    fprintf("%.4e ", iter_err(j));
end
fprintf("\n");

fprintf("num quad points in each iter\n");
for j = 1 : length(iter_err)
    fprintf("%d ", out_afom_t.num_quadpoints(j));
end
fprintf("\n");


fprintf("subspace dim in each iter\n");
for j = 1 : length(iter_err)
    fprintf("%d ", out_afom_t.dim(j));
end
fprintf("\n");
