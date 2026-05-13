function [ w,H,h,breakdown,accuracy_flag, num_oracle ] = hm_sarnoldi_fix( S,A,m,H,start_ind,param )

accuracy_flag = 0;

global V_big;

s = size(S, 1);
SAV = zeros(s, m + 1);
H(m+1,m) = 0;
breakdown = 0;

n = size(V_big, 1);
AV_big = zeros(n, m + 1);

num_oracle = 0;
for k = start_ind:(m + 1)

    w = V_big(:,k);
    if isnumeric(A)
        w = A*w;
    else
        w = A(w);
    end
    AV_big(:, k) = w;
    SAV(:, k) = S * w;
    num_oracle = num_oracle + 1;

    c = SAV(:, 1 : (k - 1))' * SAV(:, k);
    SAV(:, k) = SAV(:, k) - SAV(:, 1 : (k - 1)) * c;
    rc = SAV(:, 1 : (k - 1))' * SAV(:, k);
    SAV(:, k) = SAV(:, k) - SAV(:, 1 : (k - 1)) * rc;
    c = c + rc;
    V_big(:, k) = V_big(:, k) - V_big(:, 1 : (k - 1)) * c;
    nrm = norm(SAV(:, k));
    SAV(:, k) = SAV(:, k) / nrm;
    V_big(:, k) = V_big(:, k) / nrm;
    if k > 1
        H(1 : (k - 1), k - 1) = c;
        H(k, k - 1) = nrm;
    end

    if k <= m
        w = w - V_big(:, 1 : k) * (H(1 : k, 1 : (k - 1)) * c);
        w = w / nrm;
        V_big(:, k + 1) = w;
    end
end

% V = V_big(:,1: (m + 1));
% norm(A*V(:,1:m) - V*H)
% norm((S*A*V)'*(S*A*V) - eye(m+1))

w = V_big(:, m + 1);
h = H(m+1,m);
H = H(1:m,1:m);

end