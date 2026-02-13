function [ m,w,H,h,breakdown,accuracy_flag ] = tarnoldi_adaptive( A, m_max, param)

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

H = zeros(m_max, m_max);

cond_tol = param.cond_tol;

global V_big;
trunc = param.truncation_length;
reo = param.reorth_number;
breakdown = 0;

n = size(V_big, 1);
s0 = 30;
S = sketching_mat(s0, n, param.sketching_mat_type);
s = s0;

SV_big = zeros(s, m_max);
SV_big(:, 1) = S * V_big(: ,1);
for j = 1 : m_max

    w = V_big(:,j);
    if isnumeric(A)
        w = A*w;
    else
        w = A(w);
    end

    i_start = max([1,j-trunc+1]);
    for r = 0:reo
        for i = i_start:j
            ip = param.inner_product(w,V_big(:,i));
            H(i,j) = H(i,j) + ip(1);
            w = w - V_big(:,i)*ip(1);
        end
    end

    H(j+1,j) = norm(w);

    if abs(H(j+1,j)) < j*eps*norm(H(1:j+1,j))
        breakdown = j;
        break
    end

    w = (1/H(j+1,j))*w;
    Sw = S * w;
    if j < m_max
        V_big(:, j + 1) = w;
        SV_big(:, j + 1) = Sw;

        if s < param.ada_sketching_size_control * (j + 1)
            S_incr = sketching_mat(s0, n, param.sketching_mat_type);
            S = [S; S_incr];
            SV_big = [SV_big;
                S_incr * V_big(:, 1 : (j + 1)), zeros(s0, m_max - (j + 1))];
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
Sw = SV_big(:, m + 1);
h = H(m+1,m);
H = H(1:m,1:m);
if ~isempty(param.last_update)
    switch param.last_update
        case "orth"
            [w,H,h] = arnoldi_last_orth_update(m, w, H, h);
            [w,H,h] = arnoldi_last_orth_update(m, w, H, h);
        case "sorth"
            [w,H,h] = arnoldi_last_sorth_update(m, w, H, h, SV_big, Sw);
            [w,H,h] = arnoldi_last_sorth_update(m, w, H, h, SV_big, S * w);
    end
end
