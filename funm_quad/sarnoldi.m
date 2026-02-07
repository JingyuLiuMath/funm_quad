function [ w,H,h,breakdown,accuracy_flag] = sarnoldi( A,m,H,s,param )

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

global S;
global V_big;
global SV_big;
global SAV_big;
H(m+1,m) = 0;
trunc = param.truncation_length;
reo = param.reorth_number;
breakdown = 0;

for k = s:m,
    
    w = V_big(:,k);
    if isnumeric(A),
        w = A*w;
    else
        w = A(w);
    end
    Sw = S * w;
    SAV_big(:, k) = Sw;
    
    sj = max([1,k-trunc+1]);
    if(k==s), sj = 1; end
    for r = 0:reo,
        for j = sj:k,
            ip = param.inner_product(Sw,SV_big(:,j));
            H(j,k) = H(j,k) + ip(1);
            Sw = Sw - SV_big(:,j) * ip(1);
            w = w - V_big(:, j) * ip(1);
        end
    end
    
    H(k+1,k) = sqrt(param.inner_product(Sw, Sw));
    
    if abs(H(k+1,k)) < k*eps*norm(H(1:k+1,k))
        breakdown = k;
        break
    end
    
    Sw = Sw / H(k + 1, k);
    w = w / H(k + 1, k);
    if k < m
        SV_big(:, k + 1) = Sw;
        V_big(:, k + 1) = w;
    end
    if param.max_restarts == 1 && (~mod(k,10) && k >= 20),
        if isa(fm,'function_handle')
            c = fm(H(1:k,1:k))*eye(k,1);
        else
            [WW,DD] = eig(H(1:k,1:k));
            ee = unit(1,k);
            c = zeros(size(ee));
            for j = 1:k
                active_nodes = diag(DD);
                subdiag = diag(H(1:k+1,1:k),-1);
                fun = @(t) param.function(DD(j,j),t) .* evalnodal(t, active_nodes, subdiag);
                I = myintegral(fun,-inf,0,'AbsTol',tol,'RelTol',tol);
                c(j) = I;
            end
            c = (WW*spdiags(c,0,k,k)/WW)*ee;
        end

        if norm(c(k-9:k)) < norm(c)*param.stopping_accuracy/2,
            accuracy_flag = 1;
            breakdown = k;
            m = k;
            break
        end
    end
    
end

h = H(m+1,m);
H = H(1:m,1:m);
