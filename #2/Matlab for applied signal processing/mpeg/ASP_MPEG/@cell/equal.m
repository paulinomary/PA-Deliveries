function r = equal(a, b)

% equal - compares to entities.  
% 
%   This is the specialized version for cell arrays.  The result is 1 if the two
%   cell arrays are exactly alike.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).


if(nargin ~= 2)
    error('equal requires two arguments');
elseif(~isa(a, 'cell') | ~isa(b, 'cell'))
    error('specialized version of equal can only handle cell array arguments');
elseif(any(size(a) ~= size(b)))
    error('a and b must have same size');
end

for i = 1 : prod(size(a)),
    if(~equal(a{i}, b{i}))
        r = false;
        return
    end
end

r = true;