function m = izigzag(v, s)

% izigzag - Inverse of zigzag.
%
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

if(~isvector(v))
    error('first argument must be a vector');
elseif(~issize(s))
    error('second argument must be a size row vector');
end

[i, j] = zigzagindices(s);

for n = 1 : length(i),
    m{i(n), j(n)} = v{n};
end
