function orth_err = check_last_sorth(V, w, S)

orth_err = norm((S * V)' * (S * w)) / norm(S * w);

end