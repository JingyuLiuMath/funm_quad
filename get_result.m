function result = get_result(f_ex, method_name, param, f, out, t)

result = struct();

result.f_ex = f_ex;
result.method = method_name;
result.param = param;

result.f = f;
result.out = out;
result.time = t;

num_it = size(out.appr, 2);
rel_err = norm(f_ex - f) / norm(f_ex);

result.num_it = num_it;
result.rel_err = rel_err;

end