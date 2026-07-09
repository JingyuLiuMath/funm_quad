function [lastv, H, h, beta, breakdown, num_oracle] = shmarnoldi_fix(S, A, m, H, start_ind, param)

global V_big;
H(m + 1, m) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
breakdown = 0;
num_oracle = 0;

W_big = zeros(size(V_big));
P_big = zeros(size(S, 1), size(V_big, 2));
Q_big = zeros(size(S, 1), size(V_big, 2));
U_big = zeros(size(S, 1), size(V_big, 2));

if isnumeric(A)
    W_big(:, 1 : start_ind) = A * V_big(:, 1 : start_ind);
else
    W_big(:, 1 : start_ind) = A(V_big(:, 1 : start_ind));
end
num_oracle = num_oracle + start_ind;

P_big(:, 1 : start_ind) = S * V_big(:, 1 : start_ind);
Q_big(:, 1 : start_ind) = S * W_big(:, 1 : start_ind);

lastv = V_big(:, start_ind);
lastp = P_big(:, start_ind);
lastw = W_big(:, start_ind);
lastq = Q_big(:, start_ind);

delta = lastq' * lastp;
if abs(delta) < eps * norm(lastq) * norm(lastp)
    breakdown = start_ind - 1;
    beta = 0;
    h = 0;
    H = H(1 : m, 1 : m);
    return
end

beta = sqrt(abs(delta));
lastv = lastv / beta;
lastp = lastp / beta;
lastw = lastw / beta;
lastq = lastq / beta;

V_big(:, start_ind) = lastv;
P_big(:, start_ind) = lastp;
W_big(:, start_ind) = lastw;
Q_big(:, start_ind) = lastq;

if start_ind > 1
    cross = Q_big(:, 1 : start_ind)' * P_big(:, 1 : start_ind);
    offdiag = cross - diag(diag(cross));
    if norm(offdiag, "fro") > 1e-10 * max(1, norm(cross, "fro"))
        error("FUNM_QUAD:shmarnoldi_fix", ...
            "shmarnoldi_fix requires the carried basis to satisfy (S*A*V)'*(S*V) diagonal.");
    end
end

for i = 1 : start_ind
    rho = Q_big(:, i)' * P_big(:, i);
    if abs(rho) < eps * norm(Q_big(:, i)) * norm(P_big(:, i))
        error("FUNM_QUAD:shmarnoldi_fix", ...
            "Biorthogonal breakdown while initializing column %d.", i);
    end
    U_big(:, i) = Q_big(:, i) / conj(rho);
end

for j = start_ind : m
    lastv = lastw;
    lastp = lastq;

    j_trunc_start = max([1, j - trunc + 1]);
    if j == start_ind
        j_trunc_start = 1;
    end
    for r = 0 : reo
        for i = j_trunc_start : j
            ip = U_big(:, i)' * lastp;
            H(i, j) = H(i, j) + ip;
            lastv = lastv - V_big(:, i) * ip;
            lastp = lastp - P_big(:, i) * ip;
        end
    end

    if isnumeric(A)
        lastw = A * lastv;
    else
        lastw = A(lastv);
    end
    lastq = S * lastw;
    num_oracle = num_oracle + 1;

    delta = lastq' * lastp;
    H(j + 1, j) = sqrt(abs(delta));
    if abs(H(j + 1, j)) < j * eps * norm(H(1 : (j + 1), j))
        breakdown = j;
        break
    end

    lastv = lastv / H(j + 1, j);
    lastp = lastp / H(j + 1, j);
    lastw = lastw / H(j + 1, j);
    lastq = lastq / H(j + 1, j);

    if j < (m + 1)
        V_big(:, j + 1) = lastv;
        P_big(:, j + 1) = lastp;
        W_big(:, j + 1) = lastw;
        Q_big(:, j + 1) = lastq;

        rho = lastq' * lastp;
        if abs(rho) < eps * norm(lastq) * norm(lastp)
            error("FUNM_QUAD:shmarnoldi_fix", ...
                "Biorthogonal breakdown at column %d.", j + 1);
        end
        U_big(:, j + 1) = lastq / conj(rho);
    end
end

lastv = V_big(:, m + 1);
V_big(:, m + 1) = 0;
h = H(m + 1, m);
H = H(1 : m, 1 : m);

end