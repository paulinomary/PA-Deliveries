function r = equal(a, b)

% equal - compares to entities.  
% 
%   This is the specialized version for int32 arrays.  The result is 1 if the 
%   two int32 arrays are exactly alike.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).


if(nargin ~= 2)
    error('equal requires two arguments');
elseif(~isa(a, 'int32') | ~isa(b, 'int32'))
    error('specialized version of equal can only handle int32 array arguments');
elseif(any(size(a) ~= size(b)))
    error('a and b must have same size');
end

r = all(a(:) == b(:));
