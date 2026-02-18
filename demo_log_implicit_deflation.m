% Use the quadrature-based restart algorithm FUNM_QUAD described in
%
%  A. Frommer, S. G\"{u}ttel, and M. Schweitzer: Efficient and
%  stable Arnoldi restarts for matrix functions based on quadrature,
%  SIAM J. Matrix Anal. Appl., 35:661--683, 2014.
%
% to compute log(I+A)/A*b for a 40x40 finite difference discretization of
% the two-dimensional Laplace operator both with and without implicit
% deflation.
%%
clear;
close all;
rng(2026);

quad_tol = 1e-7;
stop_tol = 1e-8;

m = 30;
max_restarts = 15;
ada_max_restarts = 30;
truncation_length = 0;
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
        method = "afomt";
    case 3
        method = "asfomt";
end
f = @(x) log(1+x)/x;
N = 40;
e = ones(N,1);
A = (N+1)^2*gallery('poisson',N);
b = kron(e,e);
b = b/norm(b);

load exact_solutions

%% Choose parameters for the FUNM_QUAD algorithm (no implicit deflation)
addpath('funm_quad')
param.function = 'log';
param.restart_length = m;              % Each restart cycle consists of 20 Arnoldi iterations
param.max_restarts = max_restarts;                % We perform at most 20 restart cycles
param.tol = quad_tol;                      % Tolerance for quadrature rule
param.hermitian = 0;                    % set 0 if A is not Hermitian
param.V_full = 0;                       % set 1 if you need Krylov basis
param.H_full = 0;                       % do not store all Hessenberg matrices
param.exact = exact_poisson_log;        % Exact solution. If not known set to []
param.stopping_accuracy = stop_tol;        % stopping accuracy
param.inner_product = @(a,b) b'*a;      % Use standard euclidean inner product
param.thick = [];                       % No implicit deflation is performed
param.min_decay = 0.95;                 % we desire linear error reduction of rate < .95
param.waitbar = 0;                      % show waitbar
param.reorth_number = 0;                % #reorthogonalizations
param.truncation_length = inf;          % truncation length for Arnoldi
param.verbose = 2;                      % print information about progress of algorithm

%% compute log(I+A)/A*b using FUNM_QUAD without implicit deflation
if method == "benchmark"
    tic
    [f1,out1] = funm_quad(A,b,param);
    toc
elseif method == "afomt"
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
    toc
elseif method == "asfomt"
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
    toc
end


% plot convergence curve and number of quadrature points
semilogy(out1.err,'g--+')
hold on
for k = 2:length(out1.err)
    text(k+0.1,2*out1.err(k),num2str(out1.num_quadpoints(k)),'Color',[0 1 0],'FontSize',16,'Rotation',45);
end

%% adapt parameters for FUNM_QUAD algorithm (with with implicit deflation)
param.thick = @thick_quad;              % Thick restart function for implicit deflation
param.number_thick = 5;                 % Number of target eigenvalues for implicit deflation

if method == "benchmark"
    tic
    [f2,out2] = funm_quad(A,b,param);
    toc
elseif method == "afomt"
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
    toc
elseif method == "asfomt"
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
    toc
end

% plot convergence curve and number of quadrature points
semilogy(out2.err,'m--+')
for k = 2:length(out2.err)
    text(k+0.1,2*out2.err(k),num2str(out2.num_quadpoints(k)),'Color',[1 0 1],'FontSize',16,'Rotation',45);
end

legend('without deflation','implicit deflation')
xlabel('cycle'); ylabel('absolute 2-norm error')
hold off

max_iter = max(size(out1.appr, 2), size(out2.appr, 2));
xticks(1 : ceil(max_iter / 10) : max_iter);
if method == "benchmark"
    file_name = "implicit_deflation_rel_err_" + method + "_" + string(N) + "_" + string(m);
    saveas(gcf, "./figure/implicit_deflation/" + file_name + ".eps", "epsc");
else
    file_name = "implicit_deflation_rel_err_" + method + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/implicit_deflation/" + file_name + ".eps", "epsc");
    
    figure();
    plot(out1.dim, '--o', "DisplayName", "without deflation");
    hold on;
    plot(out2.dim, '--p', "DisplayName", "implicit deflation");
    legend;
    xticks(1 : ceil(max_iter / 10) : max_iter);
    yticks(unique([out1.dim, out2.dim]));
    xlabel('cycle');
    ylabel('subspace dim');
    file_name = "implicit_deflation_subspace_dim_" + method + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/implicit_deflation/" + file_name + ".eps", "epsc");
end