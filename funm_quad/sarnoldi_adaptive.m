function [ m,w,H,h,breakdown,accuracy_flag,check_result ] = sarnoldi_adaptive( A, m_max, H, start_ind, param)

accuracy_flag = 0;
check_result = struct();
check_result.before = struct('rel_err_AD', [], 'rel_orth_err', [], 'cond_V', []);
check_result.after = struct('rel_err_AD', [], 'rel_orth_err', [], 'cond_V', []);
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
global S
H(m_max + 1, m_max) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
cond_tol = param.cond_tol;
breakdown = 0;

if isempty(S)
    n = size(V_big, 1);
    if param.sketching_mat_type == "exact"
        s0 = n;
    else
        s0 = 30;
    end
    S = sketching_mat(s0, n, param.sketching_mat_type);
end
s = size(S, 1);

SV_big = zeros(s, m_max);
Sw = S * V_big(:, 1);
beta = norm(Sw);
V_big(:, 1) = V_big(:, 1) / beta;
SV_big(:, 1) = Sw / beta;
for j = start_ind : m_max

    w = V_big(:,j);
    if isnumeric(A)
        w = A*w;
    else
        w = A(w);
    end
    Sw = S * w;

    i_start = max([1,j-trunc+1]);
    for r = 0:reo
        for i = i_start:j
            ip = param.inner_product(Sw,SV_big(:,i));
            H(i,j) = H(i,j) + ip(1);
            w = w - V_big(:,i)*ip(1);
            Sw = Sw - SV_big(:,i)*ip(1);
        end
    end

    H(j+1,j) = norm(Sw);

    if abs(H(j+1,j)) < j*eps*norm(H(1:j+1,j))
        breakdown = j;
        break
    end

    w = w / H(j+1,j);
    Sw = Sw / H(j+1,j);
    if j < m_max
        V_big(:, j + 1) = w;
        SV_big(:, j + 1) = Sw;

        % Sketching matrix S is fixed (global)

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
h = H(m+1,m);
H = H(1:m,1:m);
if ~isempty(param.update)
    switch param.update
        case "last_orth"
            if param.check == 1
                check_result.before.rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                check_result.before.rel_orth_err = check_last_orth(V_big(:, 1 : m), w);
                check_result.before.cond_V = check_cond_V(V_big(:, 1 : m));
            end

            % [w, H, h] = arnoldi_last_orth_update(m, w, H, h);
            [w, H, h] = arnoldi_last_orth_update(m, w, H, h);

            if param.check == 1
                check_result.after.rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                check_result.after.rel_orth_err = check_last_orth(V_big(:, 1 : m), w);
                check_result.after.cond_V = check_cond_V(V_big(:, 1 : m));
            end
        case "last_sorth"
            if param.check == 1
                check_result.before.rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                check_result.before.rel_orth_err = check_last_sorth(V_big(:, 1 : m), w, S);
                check_result.before.cond_V = check_cond_V(V_big(:, 1 : m));
            end

            % [w, H, h, Sw] = arnoldi_last_sorth_update(m, w, H, h, SV_big, Sw);
            [w, H, h, ~] = arnoldi_last_sorth_update(m, w, H, h, SV_big, Sw);

            if param.check == 1
                check_result.after.rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                check_result.after.rel_orth_err = check_last_sorth(V_big(:, 1 : m), w, S);
                check_result.after.cond_V = check_cond_V(V_big(:, 1 : m));
            end
        case "whitening"
            if param.check == 1
                check_result.before.rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                check_result.before.rel_orth_err = check_last_sorth(V_big(:, 1 : m), w, S);
                check_result.before.cond_V = check_cond_V(V_big(:, 1 : m));
            end

            [w, H, h, ~, ~] = arnoldi_whitening_update(m, w, H, h, SV_big, Sw);

            if param.check == 1
                check_result.after.rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
                check_result.after.rel_orth_err = check_last_sorth(V_big(:, 1 : m), w, S);
                check_result.after.cond_V = check_cond_V(V_big(:, 1 : m));
            end
    end
else
    if param.check == 1
        check_result.before.rel_err_AD = check_arnoldi(m, V_big(:, 1 : m), w, H, h, A);
        check_result.before.rel_orth_err = check_last_orth(V_big(:, 1 : m), w);
        check_result.before.cond_V = check_cond_V(V_big(:, 1 : m));
    end
end

V_big(:, (m + 1) : end) = 0;
