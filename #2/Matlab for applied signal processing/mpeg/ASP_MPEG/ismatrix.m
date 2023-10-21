function r = ismatrix(m)

% ismatrix - True for matrices.
%
%   r = ismatrix(m)  Returns true if 'm' is a matrix, that is, an array with
%   1 or 2 dimensions.
%   Scalars and vectors are also matrices.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

r = ndims(m) == 2;