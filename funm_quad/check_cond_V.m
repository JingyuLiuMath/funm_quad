function check_cond_V(V)

cond_V = cond(V);

if cond_V > 1e6
    fprintf("cond(V): %.4e\n", cond(V));
    % keyboard;
end

end