function [lastv, H, h, beta, breakdown, num_oracle] = arnoldi_fix(A, m, H, start_ind, param)

global V_big;
H(m + 1, m) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
breakdown = 0;
num_oracle = 0;

lastv = V_big(:, start_ind);
beta = norm(lastv);
lastv = lastv / beta;
V_big(:, start_ind) = lastv;

for j = start_ind : m
    if isnumeric(A)
        lastv = A * lastv;
    else
        lastv = A(lastv);
    end
    num_oracle = num_oracle + 1;

    j_trunc_start = max([1, j - trunc + 1]);
    if j == start_ind
        j_trunc_start = 1;
    end
    for r = 0 : reo
        for i = j_trunc_start : j
            ip = V_big(:, i)' * lastv;
            H(i, j) = H(i, j) + ip;
            lastv = lastv - V_big(:, i) * ip;
        end
    end

    H(j + 1, j) = norm(lastv);
    if abs(H(j + 1, j)) < j * eps * norm(H(1 : (j + 1), j))
        breakdown = j;
        break
    end

    lastv = lastv / H(j + 1, j);
    if j < m
        V_big(:, j + 1) = lastv;
    end
end

h = H(m + 1, m);
H = H(1 : m, 1 : m);


end