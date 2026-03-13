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

method_list = ["benchmark", ...
    "fix FOM-t whitening", ...
    "fix FOM-s", ...
    "ada FOM-t last-orth", "ada FOM-t last-sorth", "ada FOM-t whitening", ...
    "ada FOM-s whitening", "ada FOM-st whitening"];
method_name = "ada FOM-t last-sorth";

f = @(x) log(1+x)/x;
N = 40;
e = ones(N,1);
A = (N+1)^2*gallery('poisson',N);
b = kron(e,e);
b = b/norm(b);

load exact_solutions
f_ex = exact_poisson_log;

%% Choose parameters for the FUNM_QUAD algorithm (no implicit deflation)
addpath('funm_quad')
basic_param.function = 'log';
basic_param.restart_length = m;              % Each restart cycle consists of 20 Arnoldi iterations
basic_param.max_restarts = max_restarts;                % We perform at most 20 restart cycles
basic_param.tol = quad_tol;                      % Tolerance for quadrature rule
basic_param.hermitian = 0;                    % set 0 if A is not Hermitian
basic_param.V_full = 0;                       % set 1 if you need Krylov basis
basic_param.H_full = 0;                       % do not store all Hessenberg matrices
basic_param.exact = [];        % Exact solution. If not known set to []
basic_param.stopping_accuracy = stop_tol;        % stopping accuracy
basic_param.inner_product = @(a,b) b'*a;      % Use standard euclidean inner product
basic_param.thick = [];                       % No implicit deflation is performed
basic_param.min_decay = 0.95;                 % we desire linear error reduction of rate < .95
basic_param.waitbar = 0;                      % show waitbar
basic_param.reorth_number = 0;                % #reorthogonalizations
basic_param.truncation_length = inf;          % truncation length for Arnoldi
basic_param.verbose = 2;                      % print information about progress of algorithm

add_param = construct_ada_param(...
    truncation_length, ...
    max_num_quad_points, ...
    sketching_mat_type, sketching_size, ...
    ada_max_restarts, ada_sketching_size_control, cond_tol);

%% compute log(I+A)/A*b using FUNM_QUAD without implicit deflation
result = run_single_method(A, b, f_ex, method_name, basic_param, add_param);

f1 = result.f;
out1 = result.out;

%% adapt parameters for FUNM_QUAD algorithm (with with implicit deflation)
basic_param.thick = @thick_quad;              % Thick restart function for implicit deflation
basic_param.number_thick = 5;                 % Number of target eigenvalues for implicit deflation

result = run_single_method(A, b, f_ex, method_name, basic_param, add_param);
f2 = result.f;
out2 = result.out;

%% plot
figure();
err1 = vecnorm(f_ex - out1.appr) / norm(f_ex);
semilogy(err1,'g--+')
hold on
for k = 2:length(err1)
    text(k+0.1,2*err1(k),num2str(out1.num_quadpoints(k)),'Color',[0 1 0],'FontSize',16,'Rotation',45);
end
err2 = vecnorm(f_ex - out2.appr) / norm(f_ex);
semilogy(err2,'m--+')
for k = 2:length(err2)
    text(k+0.1,2*err2(k),num2str(out2.num_quadpoints(k)),'Color',[1 0 1],'FontSize',16,'Rotation',45);
end

legend('without deflation','implicit deflation')
xlabel('cycle'); ylabel('absolute 2-norm error')
hold off

max_iter = max(size(out1.appr, 2), size(out2.appr, 2));
xticks(1 : ceil(max_iter / 10) : max_iter);
prefix_char = char(method_name);
ada_flag = strcmp(prefix_char(1:3), 'ada');
if ada_flag == 0
    file_name = "implicit_deflation_rel_err_" + method_name + "_" + string(N) + "_" + string(m);
    saveas(gcf, "./figure/implicit_deflation/" + file_name + ".eps", "epsc");
else
    file_name = "implicit_deflation_rel_err_" + method_name + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
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
    file_name = "implicit_deflation_subspace_dim_" + method_name + "_" + string(N) + "_" + string(m) + "_" + string(truncation_length);
    saveas(gcf, "./figure/implicit_deflation/" + file_name + ".eps", "epsc");
end