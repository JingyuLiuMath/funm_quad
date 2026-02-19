function [w, Sw, H, h] = arnoldi_last_sorth_update(m, w, H, h, SV_big, Sw)

global V_big;

c = SV_big \ Sw;
w = w - V_big(:, 1 : m) * c;
Sw = Sw - SV_big * c;
H(1:m, m) = H(1:m, m) + c * h;
norm_w = norm(w);
w = w / norm_w;
h = h * norm_w;

end