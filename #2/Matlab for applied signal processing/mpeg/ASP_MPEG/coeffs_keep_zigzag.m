function dr = coeffs_keep_zigzag(d,N)

v = zigzag8(d);
ve = [v(1:N),zeros(1,64-N)];
dr = izigzag8(ve);
