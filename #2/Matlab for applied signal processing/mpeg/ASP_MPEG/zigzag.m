function v = zigzag(m)

% zigzag - Returns a vector with the elements of a matrix obtained in zigzag
% order (as the DCT coeficients in JPEG).
%
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

if(~ismatrix(m))
    error('argument must be a matrix');
end

[i, j] = zigzagindices(size(m));

for n = 1 : length(i),
    v(n) = m(i(n), j(n));
end
