function [m, lastv, H, h, beta, breakdown, num_oracle, S] = arnoldi_adaptive(S, A, m_max, H, start_ind, param)

global V_big;
H(m_max + 1, m_max) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
breakdown = 0;
num_oracle = 0;

n = size(V_big, 1);
s0 = 30;
if isempty(S)
    S = sketching_mat(s0, n, param.sketching_mat_type);
end

Q_big = zeros(size(S, 1), size(V_big, 2));

lastv = V_big(:, start_ind);
beta = norm(lastv);
lastv = lastv / beta;
V_big(:, start_ind) = lastv;

Q_big(:, 1 : start_ind) = S * V_big(:, 1 : start_ind);
lastq = Q_big(:, start_ind);

examine_period = min(ceil((m_max - start_ind + 1) / 10), 5);
examine_period = max(examine_period, 1);

for j = start_ind : m_max
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
    lastq = S * lastv;

    if j < m_max
        V_big(:, j + 1) = lastv;
        Q_big(:, j + 1) = lastq;

        if size(S, 1) < param.sketching_size_control * (j + 1)
            S_incr = sketching_mat(s0, n, param.sketching_mat_type);
            S = [S; S_incr];
            Q_big = [Q_big; S_incr * V_big(:, 1 : (j + 1)), zeros(size(S_incr, 1), size(Q_big, 2) - j - 1)];
            lastq = Q_big(:, j + 1);
            if param.verbose >= 2
                fprintf("sketching size increased to %d\n", size(S, 1));
            end
        end

        if mod(num_oracle, examine_period) == 0
            if cond(Q_big(:, 1 : (j + 1))) > param.cond_tol
                break
            end
        end
    end
end

m = j;
h = H(m + 1, m);
H = H(1 : m, 1 : m);
SV_big = Q_big(:, 1 : m);
V_big(:, m + 1) = 0;

if ~isempty(param.update)
    switch param.update
        case "last_orth"
            print_check_metrics("before update", m, lastv, H, h, A, V_big(:, 1 : m), [], param);
            [lastv, H, h] = arnoldi_last_orth_update(m, lastv, H, h);
            print_check_metrics("after update", m, lastv, H, h, A, V_big(:, 1 : m), [], param);
        case "last_sorth"
            print_check_metrics("before update", m, lastv, H, h, A, V_big(:, 1 : m), S, param);
            [lastv, H, h, ~] = arnoldi_last_sorth_update(m, lastv, H, h, SV_big, lastq);
            print_check_metrics("after update", m, lastv, H, h, A, V_big(:, 1 : m), S, param);
        case "last_hmorth"
            if should_check(param)
                AV_big = compute_AV(A, V_big(:, 1 : m));
                print_check_metrics("before update", m, lastv, H, h, A, AV_big, [], param);
            end
            [lastv, H, h] = arnoldi_last_hmorth_update(m, lastv, H, h);
            if should_check(param)
                print_check_metrics("after update", m, lastv, H, h, A, AV_big, [], param);
            end
        case "last_shmorth"
            if should_check(param)
                AV_big = compute_AV(A, V_big(:, 1 : m));
                print_check_metrics("before update", m, lastv, H, h, A, AV_big, S, param);
            end
            [lastv, H, h] = arnoldi_last_shmorth_update(m, lastv, H, h, SV_big, lastq);
            if should_check(param)
                print_check_metrics("after update", m, lastv, H, h, A, AV_big, S, param);
            end
        otherwise
            error("FUNM_QUAD: unknown adaptive Arnoldi update.");
    end
else
    print_check_metrics("no update", m, lastv, H, h, A, V_big(:, 1 : m), [], param);
end

end

function print_check_metrics(label, m, lastv, H, h, A, orth_basis, S, param)

if should_check(param)
    global V_big;
    
    rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), lastv, H, h, A);
    if isempty(S)
        rel_orth_err = check_last_orth(orth_basis, lastv);
    else
        rel_orth_err = check_last_sorth(orth_basis, lastv, S);
    end
    cond_V = check_cond_V(V_big(:, 1 : m));

    fprintf("%s:\n", label);
    fprintf("  rel_err_AD: %.1e\n", rel_err_AD);
    fprintf("  rel_orth_err: %.1e\n", rel_orth_err);
    fprintf("  cond_V: %.1e\n", cond_V);
end

end

function flag = should_check(param)

flag = isfield(param, "check") && param.check == 1;

end

function AV_big = compute_AV(A, V)

if isnumeric(A)
    AV_big = A * V;
else
    AV_big = A(V);
end

end