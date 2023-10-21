function r = equal(a, b)

% equal - compares to entities.  
% 
%   This is the specialized version for single arrays.  The result is 1 if the 
%   two double arrays are exactly alike.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).


if(nargin ~= 2)
    error('equal requires two arguments');
elseif(~isa(a, 'single') | ~isa(b, 'single'))
    error('specialized version of equal can only handle single array arguments');
elseif(any(size(a) ~= size(b)))
    error('a and b must have same size');
end

r = all(a(:) == b(:));
