function m = izigzag8(v)

% izigzag8 - Inverse of zigzag8.
%
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

if(~isvector(v) | length(v) ~= 64)
    error('argument must be a vector with 64 elements');
end

i = [1 1 2 3 2 1 1 2 3 4 5 4 3 2 1 1 2 3 4 5 6 7 6 5 4 3 2 1 ... 
     1 2 3 4 5 6 7 8 8 7 6 5 4 3 2 3 4 5 6 7 8 8 7 6 5 4 5 6 7 8 8 7 6 7 8 8];

j = [1 2 1 1 2 3 4 3 2 1 1 2 3 4 5 6 5 4 3 2 1 1 2 3 4 5 6 7 ...
     8 7 6 5 4 3 2 1 2 3 4 5 6 7 8 8 7 6 5 4 3 4 5 6 7 8 8 7 6 5 6 7 8 8 7 8];

for n = 1 : length(i),
    m(i(n), j(n)) = v(n);
end
