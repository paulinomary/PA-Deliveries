function dr = coeffs_keep_higher(d,N)

d2 = d(:);
[v,i] = sort(abs(d2),1,'descend');
v = d2(i); ii(i) = 1:64;
ve = [v(1:N);zeros(64-N,1)];
veu = ve(ii);
dr = reshape(veu,8,8);
