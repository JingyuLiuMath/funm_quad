function [lastv, H, h, beta, breakdown, num_oracle] = shmarnoldi_fix(S, A, m, H, start_ind, param)

global V_big;
H(m + 1, m) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
breakdown = 0;
num_oracle = 0;

P_big = zeros(size(S, 1), size(V_big, 2));
Q_big = zeros(size(S, 1), size(V_big, 2));
U_big = zeros(size(S, 1), size(V_big, 2));

P_big(:, 1 : start_ind) = S * V_big(:, 1 : start_ind);

if start_ind > 1
    if isnumeric(A)
        AV = A * V_big(:, 1 : (start_ind - 1));
    else
        AV = A(V_big(:, 1 : (start_ind - 1)));
    end
    num_oracle = num_oracle + start_ind - 1;
    Q_big(:, 1 : (start_ind - 1)) = S * AV;
end

lastv = V_big(:, start_ind);
lastp = P_big(:, start_ind);

beta = norm(lastp);
if beta < eps * norm(lastv)
    breakdown = start_ind - 1;
    beta = 0;
    h = 0;
    H = H(1 : m, 1 : m);
    return
end

lastv = lastv / beta;
lastp = lastp / beta;

V_big(:, start_ind) = lastv;
P_big(:, start_ind) = lastp;

if start_ind > 1
    cross = Q_big(:, 1 : (start_ind - 1))' * P_big(:, 1 : start_ind);
    target = [diag(diag(cross(:, 1 : (start_ind - 1)))), zeros(start_ind - 1, 1)];
    if norm(cross - target, "fro") > 1e-10 * max(1, norm(cross, "fro"))
        error("FUNM_QUAD:shmarnoldi_fix", ...
            "shmarnoldi_fix requires the carried basis to satisfy (S*A*V)'*(S*V) upper orthogonality.");
    end
end

for i = 1 : (start_ind - 1)
    rho = Q_big(:, i)' * P_big(:, i);
    if abs(rho) < eps * norm(Q_big(:, i)) * norm(P_big(:, i))
        error("FUNM_QUAD:shmarnoldi_fix", ...
            "Biorthogonal breakdown while initializing column %d.", i);
    end
    U_big(:, i) = Q_big(:, i) / conj(rho);
end

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

    % Scale the residual with its sketched norm.  This avoids computing
    % S*A*V(:,j+1) at the end of the last step.
    H(j + 1, j) = norm(lastp);
    if abs(H(j + 1, j)) < j * eps * norm(H(1 : (j + 1), j))
        breakdown = j;
        break
    end

    lastv = lastv / H(j + 1, j);
    lastp = lastp / H(j + 1, j);

    if j < m
        V_big(:, j + 1) = lastv;
        P_big(:, j + 1) = lastp;
    end
end

h = H(m + 1, m);
H = H(1 : m, 1 : m);

end
