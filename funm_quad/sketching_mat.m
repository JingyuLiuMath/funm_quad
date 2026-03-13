function S = sketching_mat(s, N, sketching)

switch sketching
    case "Gaussian"
        S = randn(s, N);
    case "Clarkson-Woodruff"
        S = clarkson_woodruff(s, N);
    case "sparse sign"
        S = sparse_sign(s, N);
    case "exact"
        assert(s == N, "s must be equal to N!");
        S = speye(N);
end

end
