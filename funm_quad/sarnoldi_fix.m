function [lastv, H, h, beta, breakdown, num_oracle] = sarnoldi_fix(S, A, m, H, start_ind, param)

global V_big;
H(m + 1, m) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
breakdown = 0;
num_oracle = 0;

P_big = zeros(size(S, 1), size(V_big, 2));
P_big(:, 1 : start_ind) = S * V_big(:, 1 : start_ind);

lastp = P_big(:, start_ind);
lastv = V_big(:, start_ind);

beta = norm(lastp);
lastp = lastp / beta;
lastv = lastv / beta;

P_big(:, start_ind) = lastp;
V_big(:, start_ind) = lastv;

for j = start_ind : m
    if isnumeric(A)
        lastv = A * lastv;
    else
        lastv = A(lastv);
    end
    lastp = S * lastv;
    num_oracle = num_oracle + 1;

    j_trunc_start = max([1, j - trunc + 1]);
    for r = 0 : reo
        for i = j_trunc_start : j
            ip = P_big(:, i)' * lastp;
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