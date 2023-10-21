function r = equal(a, b)

% equal - compares to entities.  
% 
%   This is the specialized version for sparse matrices.  The result is 1 if the 
%   two sparse matrices are exactly alike.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).


if(nargin ~= 2)
    error('equal requires two arguments');
elseif(~isa(a, 'sparse') | ~isa(b, 'sparse'))
    error('specialized version of equal can only handle sparse matrix arguments');
elseif(any(size(a) ~= size(b)))
    error('a and b must have same size');
end

r = full(all(a(:) == b(:)));
