function [w, H, h, Sw] = arnoldi_last_shmorth_update(m, w, H, h, SV_big, SAV_big, Sw)

global V_big;

% SAV_big is kept in the signature for compatibility.  Using the QR
% factorization below avoids forming the sketched normal equations.
em = zeros(m,1); em(m)=1;
[~,R] = qr([SV_big,Sw], "econ");

R11 = R(1:m,1:m);
r   = R(1:m,end);
rho = R(end,end);
M = H'*R11' + h*(em*r');
c = R11\(r + M \ (h*rho^2*em));

w = w - V_big(:, 1 : m) * c;
Sw = Sw - SV_big * c;
H(1:m, m) = H(1:m, m) + c * h;
norm_w = norm(Sw);
w = w / norm_w;
h = h * norm_w;

end
