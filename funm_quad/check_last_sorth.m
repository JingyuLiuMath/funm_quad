function check_last_sorth(V, w, S)

orth_err = norm((S * V)' * (S * w)) / norm(S * w);

if orth_err > 1e-10
    fprintf("rel orth err: %.4e\n", orth_err);
    % keyboard;
end

end