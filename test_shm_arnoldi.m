rng(1);
N = 1000;
if 0
    A = randn(N);
    b = randn(N, 1);
else
    A = complex(randn(N), randn(N));
    b = complex(randn(N, 1), randn(N, 1));
end

debug_mode = 1;

m = 1;
N = size(b, 1);

% s = 2 * m + 1;
% S = randn(s, N) / sqrt(s);
s = N;
S = speye(N);

V = zeros(N, m + 1);
P = zeros(s, m + 1);
Q = zeros(s, m);
Z = zeros(s, m);
H = zeros(m + 1, m);
num_oracle = 0;

v = b;
p = S * v;
beta = norm(p);
if beta < eps * norm(b)
    error("test_lab:breakdown", ...
        "Normalization breakdown at the first vector.");
end

p = p / beta;
v = v / beta;

V(:, 1) = v;
P(:, 1) = p;

for j = 1 : m
    v = A * v;
    p = S * v;
    num_oracle = num_oracle + 1;
    Q(:, j) = p;

    rho = Q(:, j)' * P(:, j);
    if abs(rho) < eps * norm(Q(:, j)) * norm(P(:, j))
        error("test_lab:breakdown", ...
            "Biorthogonal breakdown at column %d.", j);
    end
    Z(:, j) = Q(:, j) / conj(rho);
    
    for i = 1 : j
        if debug_mode == 1
            denom = Z(:, i)' * P(:, i);
            if abs(denom - 1) > 1e-10
                error("test_lab:breakdown", ...
                    "Biorthogonal normalization lost at i = %d, j = %d.", i, j);
            end
        end

        H(i, j) = Z(:, i)' * p;
        p = p - P(:, i) * H(i, j);
        v = v - V(:, i) * H(i, j);
    end

    if debug_mode == 1
        orth_err = norm(Q(:, 1 : j)' * p);
        biorth_err = norm(Z(:, 1 : j)' * p);
        fprintf("j = %2d, orth_err before scaling: %.1e\n", j, orth_err);
        fprintf("j = %2d, biorth_err before scaling: %.1e\n", j, biorth_err);
    end

    % Scale the residual with its sketched norm. This avoids computing
    % S*A*V(:,j+1) at the end of the last step.
    H(j + 1, j) = norm(p);
    if abs(H(j + 1, j)) < j * eps * norm(H(1 : (j + 1), j))
        error("test_lab:breakdown", ...
            "Normalization breakdown at j = %d.", j);
    end

    p = p / H(j + 1, j);
    v = v / H(j + 1, j);

    P(:, j + 1) = p;
    V(:, j + 1) = v;
end

%% check
V = V(:, 1 : (m + 1));
H = H(1 : (m + 1), 1 : m);
P = P(:, 1 : (m + 1));
Q = Q(:, 1 : m);
Z = Z(:, 1 : m);

AV = A * V(:, 1 : m);
arnoldi_decomp_err = norm(AV - V * H, "fro") / norm(AV, "fro");
fprintf("arnoldi_decomp_err: %.1e\n", arnoldi_decomp_err);

SV_consistency_err = norm(P - S * V, "fro");
fprintf("SV_consistency_err: %.1e\n", SV_consistency_err);

SAV = S * AV;
Q_consistency_err = norm(Q - SAV, "fro");
fprintf("Q_consistency_err: %.1e\n", Q_consistency_err);

max_orth_err = 0;
for j = 1 : m
    max_orth_err = max(max_orth_err, norm(SAV(:, 1 : j)' * P(:, j + 1)));
end
fprintf("max_orth_err: %.1e\n", max_orth_err);

diag_err = norm(diag(Z' * P(:, 1 : m)) - ones(m, 1));
fprintf("diag_err: %.1e\n", diag_err);

fprintf("num_oracle: %d\n", num_oracle);
