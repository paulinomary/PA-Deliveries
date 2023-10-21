function r = isrow(v)

% isrow - True for row vectors.
%
%   r = isrow(v)  Returns true if 'v' is a row vector, that is, an array with
%   1 or 2 dimensions where the size of the first dimension must be one.
%   Scalars are also row vectors: isscalar(v) => isrow(v) => isvector(v).
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

r = ndims(v) == 2 & (size(v, 1) == 1);