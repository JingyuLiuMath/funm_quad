function [lastv, H, h, beta, breakdown, num_oracle] = shmarnoldi_fix(S, A, m, H, start_ind, param)

global V_big;
assert(start_ind == 1, ...
    "FUNM_QUAD:shmarnoldi_fix: shmarnoldi_fix does not support thick restart.");
H(m + 1, m) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
breakdown = 0;
num_oracle = 0;

P_big = zeros(size(S, 1), size(V_big, 2));
Q_big = zeros(size(S, 1), size(V_big, 2));
U_big = zeros(size(S, 1), size(V_big, 2));

P_big(:, 1) = S * V_big(:, 1);

lastv = V_big(:, 1);
lastp = P_big(:, 1);

beta = norm(lastp);
lastv = lastv / beta;
lastp = lastp / beta;

V_big(:, 1) = lastv;
P_big(:, 1) = lastp;

for j = start_ind : m
    if isnumeric(A)
        lastw = A * lastv;
    else
        lastw = A(lastv);
    end
    lastq = S * lastw;
    num_oracle = num_oracle + 1;
    Q_big(:, j) = lastq;

    rho = Q_big(:, j)' * P_big(:, j);
    if abs(rho) < eps * norm(Q_big(:, j)) * norm(P_big(:, j))
        error("FUNM_QUAD:shmarnoldi_fix", ...
            "Biorthogonal breakdown at column %d.", j);
    end
    U_big(:, j) = Q_big(:, j) / conj(rho);

    lastv = lastw;
    lastp = lastq;

    j_trunc_start = max([1, j - trunc + 1]);
    for r = 0 : reo
        for i = j_trunc_start : j
            ip = U_big(:, i)' * lastp;
            H(i, j) = H(i, j) + ip;
            lastp = lastp - P_big(:, i) * ip;
            lastv = lastv - V_big(:, i) * ip;
        end
    end

    H(j + 1, j) = norm(lastp);
    if abs(H(j + 1, j)) < j * eps * norm(H(1 : (j + 1), j))
        breakdown = j;
        break
    end

    lastp = lastp / H(j + 1, j);
    lastv = lastv / H(j + 1, j);
    if j < m
        P_big(:, j + 1) = lastp;
        V_big(:, j + 1) = lastv;
    end
end

h = H(m + 1, m);
H = H(1 : m, 1 : m);

end