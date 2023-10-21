function r = equal(a, b)

% equal - compares to entities.  
% 
%   This is the specialized version for uint16 arrays.  The result is 1 if the 
%   two uint16 arrays are exactly alike.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).


if(nargin ~= 2)
    error('equal requires two arguments');
elseif(~isa(a, 'uint16') | ~isa(b, 'uint16'))
    error('specialized version of equal can only handle uint16 array arguments');
elseif(any(size(a) ~= size(b)))
    error('a and b must have same size');
end

r = all(a(:) == b(:));
