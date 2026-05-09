function [w, H, h, Sw] = arnoldi_last_hm_sorth_update(m, w, H, h, SV_big, SAV_big, Sw)

global V_big;

c = (SAV_big' * SV_big) \ (SAV_big' * Sw);
w = w - V_big(:, 1 : m) * c;
Sw = Sw - SV_big * c;
H(1:m, m) = H(1:m, m) + c * h;
norm_w = norm(Sw);
w = w / norm_w;
h = h * norm_w;

end