function [ m,w,H,h,breakdown,accuracy_flag,num_oracle,S ] = arnoldi_adaptive( S, A, m_old, H, start_ind, param)

accuracy_flag = 0;
fm = 0;
tol = param.tol;
if param.max_restarts == 1
    if strcmp(param.function,'invSqrt')
        fm = @(X) inv(sqrtm(X));
    end
    if strcmp(param.function,'log')
        fm = @(X) logm(X);
    end
    if strcmp(param.function,'exp')
        fm = @(X) expm(X);
    end
end

global V_big;
H(m_old + 1, m_old) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
cond_tol = param.cond_tol;
breakdown = 0;

m_max = param.restart_length;
if isempty(S)
    s0 = 30;
    n = size(V_big, 1);
    if param.sketching_mat_type == "exact"
        S = sketching_mat(n, n, param.sketching_mat_type);
    else
        s0 = 30;
        S = sketching_mat(s0, n, param.sketching_mat_type);
    end
end
s = size(S, 1);

SV_big = S * V_big(:, 1 : start_ind);
num_oracle = 0;
for j = start_ind : m_max

    w = V_big(:,j);
    if isnumeric(A)
        w = A*w;
    else
        w = A(w);
    end
    num_oracle = num_oracle + 1;

    i_start = max([1,j-trunc+1]);
    for r = 0:reo
        for i = i_start:j
            ip = param.inner_product(w, V_big(:,i));
            H(i,j) = H(i,j) + ip(1);
            w = w - V_big(:,i)*ip(1);
        end
    end

    H(j+1,j) = norm(w);

    if abs(H(j+1,j)) < j*eps*norm(H(1:j+1,j))
        breakdown = j;
        break
    end

    w = w / H(j+1,j);
    Sw = S * w;
    if j < m_max
        V_big(:, j + 1) = w;
        SV_big(:, j + 1) = Sw;

        if s < param.sketching_size_control * (j + 1)
            S_incr = sketching_mat(s0, n, param.sketching_mat_type);
            S = [S; S_incr];
            SV_big = [SV_big;
                S_incr * V_big(:, 1 : (j + 1))];
            Sw = SV_big(:, j + 1);
            s = s + s0;
        end

        if cond(SV_big(:, 1 : (j + 1))) > cond_tol
            break;
        end
    end

    if param.max_restarts == 1 && (~mod(j,10) && j >= 20),
        if isa(fm,'function_handle')
            c = fm(H(1:j,1:j))*eye(j,1);
        else
            [WW,DD] = eig(H(1:j,1:j));
            ee = unit(1,j);
            c = zeros(size(ee));
            for i = 1:j
                active_nodes = diag(DD);
                subdiag = diag(H(1:j+1,1:j),-1);
                fun = @(t) param.function(DD(i,i),t) .* evalnodal(t, active_nodes, subdiag);
                I = myintegral(fun,-inf,0,'AbsTol',tol,'RelTol',tol);
                c(i) = I;
            end
            c = (WW*spdiags(c,0,j,j)/WW)*ee;
        end

        if norm(c(j-9:j)) < norm(c)*param.stopping_accuracy/2,
            accuracy_flag = 1;
            breakdown = j;
            m = j;
            break
        end
    end
end

m = j;
SV_big = SV_big(:, 1 : m);
h = H(m+1, m);
H = H(1 : m, 1 : m);
if ~isempty(param.update)
    switch param.update
        case "last_orth"
            if param.check == 1
                rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                rel_orth_err = check_last_orth(V_big(:, 1 : m), w);
                cond_V = check_cond_V(V_big(:, 1 : m));

                fprintf("before update:\n")
                fprintf("  rel_err_AD: %.1e\n", rel_err_AD);
                fprintf("  rel_orth_err: %.1e\n", rel_orth_err);
                fprintf("  cond_V: %.1e\n", cond_V);
            end

            % [w, H, h] = arnoldi_last_orth_update(m, w, H, h);
            [w, H, h] = arnoldi_last_orth_update(m, w, H, h);

            if param.check == 1
                rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                rel_orth_err = check_last_orth(V_big(:, 1 : m), w);
                cond_V = check_cond_V(V_big(:, 1 : m));

                fprintf("after update:\n")
                fprintf("  rel_err_AD: %.1e\n", rel_err_AD);
                fprintf("  rel_orth_err: %.1e\n", rel_orth_err);
                fprintf("  cond_V: %.1e\n", cond_V);
            end
        case "last_sorth"
            if param.check == 1
                rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                rel_orth_err = check_last_sorth(V_big(:, 1 : m), w, S);
                cond_V = check_cond_V(V_big(:, 1 : m));

                fprintf("before update:\n")
                fprintf("  rel_err_AD: %.1e\n", rel_err_AD);
                fprintf("  rel_orth_err: %.1e\n", rel_orth_err);
                fprintf("  cond_V: %.1e\n", cond_V);
            end

            % [w, H, h, Sw] = arnoldi_last_sorth_update(m, w, H, h, SV_big, Sw);
            [w, H, h, ~] = arnoldi_last_sorth_update(m, w, H, h, SV_big, Sw);

            if param.check == 1
                rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                rel_orth_err = check_last_sorth(V_big(:, 1 : m), w, S);
                cond_V = check_cond_V(V_big(:, 1 : m));

                fprintf("after update:\n")
                fprintf("  rel_err_AD: %.1e\n", rel_err_AD);
                fprintf("  rel_orth_err: %.1e\n", rel_orth_err);
                fprintf("  cond_V: %.1e\n", cond_V);
            end
        case "whitening"
            if param.check == 1
                rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                rel_orth_err = check_last_sorth(V_big(:, 1 : m), w, S);
                cond_V = check_cond_V(V_big(:, 1 : m));

                fprintf("before update:\n")
                fprintf("  rel_err_AD: %.1e\n", rel_err_AD);
                fprintf("  rel_orth_err: %.1e\n", rel_orth_err);
                fprintf("  cond_V: %.1e\n", cond_V);
            end

            [w, H, h, ~, ~] = arnoldi_whitening_update(m, w, H, h, SV_big, Sw);
            if param.check == 1
                rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                rel_orth_err = check_last_sorth(V_big(:, 1 : m), w, S);
                cond_V = check_cond_V(V_big(:, 1 : m));

                fprintf("after update:\n")
                fprintf("  rel_err_AD: %.1e\n", rel_err_AD);
                fprintf("  rel_orth_err: %.1e\n", rel_orth_err);
                fprintf("  cond_V: %.1e\n", cond_V);
            end
    end
else
    if param.check == 1
        rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
        rel_orth_err = check_last_orth(V_big(:, 1 : m), w);
        cond_V = check_cond_V(V_big(:, 1 : m));

        fprintf("no update:\n")
        fprintf("  rel_err_AD: %.1e\n", rel_err_AD);
        fprintf("  rel_orth_err: %.1e\n", rel_orth_err);
        fprintf("  cond_V: %.1e\n", cond_V);
    end
end

V_big(:, (m + 1) : end) = 0;
