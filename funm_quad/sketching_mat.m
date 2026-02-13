function S = sketching_mat(s, N, sketching)

switch sketching
    case "Gaussian"
        S = randn(s, N);
    case "Clarkson-Woodruff"
        S = clarkson_woodruff(s, N);
    case "sparse sign"
        S = sparse_sign(s, N);
end

end
