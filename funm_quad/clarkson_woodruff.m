function S = clarkson_woodruff(s, N)

i = randi(s,N,1);
j = 1:N;
el = 2*randi(2,N,1)-3;
S = sparse(i,j,el,s,N,N);

end
