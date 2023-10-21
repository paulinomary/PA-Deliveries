function [i, j] = zigzagindices(s)

% zigzagindices - Returns two vectors of indices necessary to zigzag scan 
% (as in JPEG) the elements of a matrix.
%
%
%   [i, j] = zigzagindices(s)  Returns the index vectors 'i' and 'j' (lines and
%   columns) for zigzag scaning a matriz with size 's'.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

if(~isrow(s) | ~size(s, 2) == 2)
    error('s should be the size of a two-dimensional matrix');
end

i = zeros(1, 0);
j = zeros(1, 0);
for d = 1 : sum(s) - 1,
    if(d > s(2))
        dco = d - s(2);
    else
        dco = 0;
    end
    if(d > s(1))
        dli = d - s(1);
    else
        dli = 0;
    end
    if(mod(d, 2) ~= 0)
        i = cat(2, i, d - dli : -1 : 1 + dco);
        j = cat(2, j, 1 + dli : d - dco);
    else
        i = cat(2, i, 1 + dco : d - dli);
        j = cat(2, j, d - dco : -1 : 1 + dli);
    end
end 