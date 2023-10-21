function r = equal(a, b)

% equal - compares to entities.  
% 
%   This is the specialized version for int8 arrays.  The result is 1 if the 
%   two int8 arrays are exactly alike.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).


if(nargin ~= 2)
    error('equal requires two arguments');
elseif(~isa(a, 'int8') | ~isa(b, 'int8'))
    error('specialized version of equal can only handle int8 array arguments');
elseif(any(size(a) ~= size(b)))
    error('a and b must have same size');
end

r = all(a(:) == b(:));
