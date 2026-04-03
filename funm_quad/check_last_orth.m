function orth_err = check_last_orth(V, w)

orth_err = norm(V' * w) / norm(w);

end