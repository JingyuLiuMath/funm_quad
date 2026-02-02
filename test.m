clear;
close all;
rng(2026);

%% Build discretization matrix for 2D convection-diffusion problem 
N = 500;
D2 = (N+1)^2*gallery('tridiag',N);
I = speye(N);
D2 = kron(I,D2) + kron(D2,I);
o = ones(N,1);
D1 = (N+1)/2*spdiags([-o,0*o,o],-1:1,N,N);
D1 = kron(I,D1) + kron(D1,I);
A = D2 + 0*D1;

% choose time step s = 2*1e-3
s = 2*1e-3;
A = -s*A;
A = sprandn(A);
herm_err = norm(A - A', "fro")/ norm(A, "fro");
if herm_err == 0
    fprintf("not Herm mat!\n");
end
    

% choose right-hand side as normalized vector of all ones
b = ones(N^2,1); b = b/norm(b);

%% choose parameters for the FUNM_QUAD restart algorithm
addpath('funm_quad')
param.function = 'exp';
param.restart_length = 70;          % each restart cycle consists of 70 Arnoldi iterations
param.max_restarts = 15;            % perform at most 15 restart cycles
param.tol = 1e-8;                  % tolerance for quadrature rule
param.hermitian = 0;                % the matrix A is Hermitian
param.V_full = 0;                   % set 1 if you need Krylov basis
param.H_full = 0;                   % do not store all Hessenberg matrices
param.exact = [];     % exact solution. If not known set to []
param.stopping_accuracy = 1e-9;    % stopping accuracy
param.inner_product = @(a,b) b'*a;  % use standard Euclidean inner product
param.thick = [];                   % no implicit deflation is performed
param.min_decay = .95;              % we desire linear error reduction of rate < .95 
param.waitbar = 0;                  % show waitbar 
param.reorth_number = 0;            % #reorthogonalizations
param.truncation_length = inf;      % truncation length for Arnoldi 
param.verbose = 1;                  % print information about progress of algorithm

%% compute exp(A)b by quadrature-based restart algorithm
tic
[f, out] = funm_quad(A,b,param);
toc

%% compute exp(A)b by skected quadrature-based restart algorithm
param.truncation_length = 5;
tic
[sketched_f1, sketched_out1] = sketched_funm_quad_1(A,b,param);
toc

rel_err = norm(f - sketched_f1);
fprintf("rel err: %e\n", rel_err);


%% compute exp(A)b by skected quadrature-based restart algorithm 2
param.truncation_length = inf; 
tic
[sketched_f2, sketched_out2] = sketched_funm_quad_2(A,b,param);
toc

rel_err = norm(f - sketched_f2);
fprintf("rel err: %e\n", rel_err);
