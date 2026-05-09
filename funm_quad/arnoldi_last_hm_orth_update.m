function [w,H,h] = arnoldi_last_hm_orth_update(m, w, H, h, AV_big)

global V_big;

c = (AV_big' * V_big(:, 1 : m)) \ (AV_big' * w);
w = w - V_big(:, 1 : m) * c;
H(1:m, m) = H(1:m, m) + c * h;
norm_w = norm(w);
w = w / norm_w;
h = h * norm_w;

end