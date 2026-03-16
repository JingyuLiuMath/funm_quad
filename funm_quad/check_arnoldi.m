function check_arnoldi(m, V, w, H, h, A)

AV = A * V;

diff_AD = AV - (V * H + w * h * unit(m, m)');
rel_err_AD = norm(diff_AD, "fro")/ norm(AV, "fro");

if rel_err_AD > 1e-10
    fprintf("rel decomp err: %.4e\n", rel_err_AD);
    % keyboard;
end

end