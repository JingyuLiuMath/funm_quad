function [w, H, h, SV_big, Sw] = arnoldi_whitening_update(m, w, H, h, SV_big, Sw)

global V_big;

[SV_big_ex, R_ex] = qr([SV_big, Sw], "econ");
SV_big = SV_big_ex(:, 1 : m);
Sw = SV_big_ex(:, m + 1);
R = R_ex(1 : m, 1 : m);
V_ex = [V_big(:, 1 : m), w] / R_ex;
V_big(:, 1 : m) = V_ex(:, 1 : m);
w = V_ex(:, m + 1);
H_ex = (R_ex * [H; h * unit(m, m)']) / R;
H = H_ex(1:m,1:m);
h = H_ex(m + 1, m);

end