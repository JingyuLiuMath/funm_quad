function rel_err_AD = check_arnoldi(m, V, w, H, h, A)

if isnumeric(A)
    AV = A * V;
else
    AV = A(V);
end

diff_AD = AV - (V * H + w * h * unit(m, m)');
rel_err_AD = norm(diff_AD, "fro")/ norm(AV, "fro");

end