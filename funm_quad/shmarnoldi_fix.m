function [lastv, H, h, beta, breakdown, num_oracle] = shmarnoldi_fix(S, A, m, H, start_ind, param)

global V_big;
H(m + 2, m + 1) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
breakdown = 0;
num_oracle = 0;

W_big = zeros(size(V_big));
if isnumeric(A)
    W_big(:, 1 : start_ind) = A * V_big(:, 1 : start_ind);
else
    W_big(:, 1 : start_ind) = A(V_big(:, 1 : start_ind));
end
num_oracle = num_oracle + start_ind;

Q_big = zeros(size(S, 1), size(V_big, 2));
Q_big(:, 1 : start_ind) = S * W_big(:, 1 : start_ind);

lastq  = Q_big(:, start_ind);
lastw = W_big(:, start_ind);
lastv = V_big(:, start_ind);
beta = norm(lastq);
lastq = lastq / beta;
lastw = lastw / beta;
lastv = lastv / beta;
Q_big(:, start_ind) = lastq;
W_big(:, start_ind) = lastw;
V_big(:, start_ind) = lastv;

for j = start_ind : m
    lastv = lastw;
    if isnumeric(A)
        lastw = A * lastv;
    else
        lastw = A(lastv);
    end
    lastq = S * lastw;
    num_oracle = num_oracle + 1;

    j_trunc_start = max([1, j - trunc + 1]);
    if j == start_ind
        j_trunc_start = 1;
    end
    for r = 0 : reo
        for i = j_trunc_start : j
            ip = Q_big(:, i)' * lastq;
            H(i, j) = H(i, j) + ip;
            lastq = lastq - Q_big(:, i) * ip;
            lastw = lastw - W_big(:, i) * ip;
            lastv = lastv - V_big(:, i) * ip;
        end
    end

    H(j + 1, j) = norm(lastq);
    if abs(H(j + 1, j)) < j * eps * norm(H(1 : (j + 1), j))
        breakdown = j;
        break
    end

    lastq = lastq / H(j + 1, j);
    lastw = lastw / H(j + 1, j);
    lastv = lastv / H(j + 1, j);

    if j < (m + 1)
        Q_big(:, j + 1) = lastq;
        W_big(:, j + 1) = lastw;
        V_big(:, j + 1) = lastv;
    end
end

lastv = V_big(:, m + 1);
V_big(:, m + 1) = 0;
h = H(m + 1, m);
H = H(1 : m, 1 : m);

end