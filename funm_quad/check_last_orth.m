function check_last_orth(V, w)

orth_err = norm(V' * w) / norm(w);

if orth_err > 1e-10
    fprintf("rel orth err: %.4e\n", orth_err);
    % keyboard;
end

end