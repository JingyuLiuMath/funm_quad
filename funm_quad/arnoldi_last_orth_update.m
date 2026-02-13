function [w,H,h] = Arnoldi_last_orth_update(m, w, H, h)

global V_big;

c = V_big(:, 1 : m) \ w;
w = w - V_big(:, 1 : m) * c;
H(1:m, m) = H(1:m, m) + c * h;
norm_w = norm(w);
w = w / norm_w;
h = h * norm_w;

end