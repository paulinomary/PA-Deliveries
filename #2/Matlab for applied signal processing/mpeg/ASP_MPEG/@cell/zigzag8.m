function v = zigzag8(m)

% zigzag8 - Returns a vector with the elements of a matrix obtained in zigzag
% order (as the DCT coeficients in JPEG).
%
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

if(~ismatrix(m) | ~equal(size(m), [8 8]))
    error('argument must be a 8x8 matrix');
end

i = [1 1 2 3 2 1 1 2 3 4 5 4 3 2 1 1 2 3 4 5 6 7 6 5 4 3 2 1 ... 
     1 2 3 4 5 6 7 8 8 7 6 5 4 3 2 3 4 5 6 7 8 8 7 6 5 4 5 6 7 8 8 7 6 7 8 8];

j = [1 2 1 1 2 3 4 3 2 1 1 2 3 4 5 6 5 4 3 2 1 1 2 3 4 5 6 7 ...
     8 7 6 5 4 3 2 1 2 3 4 5 6 7 8 8 7 6 5 4 3 4 5 6 7 8 8 7 6 5 6 7 8 8 7 8];

for n = 1 : length(i),
    v{n} = m{i(n), j(n)};
end