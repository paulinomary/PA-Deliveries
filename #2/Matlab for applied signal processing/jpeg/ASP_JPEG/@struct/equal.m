function r = equal(a, b)

% equal - compares to entities.  
% 
%   This is the specialized version for struct arrays.  The result is 1 if the 
%   two struct arrays are exactly alike.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).


if(nargin ~= 2)
    error('equal requires two arguments');
elseif(~isa(a, 'struct') | ~isa(b, 'struct'))
    error('specialized version of equal can only handle struct array arguments');
elseif(~equal(fieldnames(a), fieldnames(b)))
    error('a and b must have same fields!');
elseif(any(size(a) ~= size(b)))
    error('a and b must have same size');
end

fields = fieldnames(a);

for i = 1 : prod(size(a)),
    for j = 1 : length(fields),
        if(~equal(getfield(a(i), fields{j}), getfield(b(i), fields{j})))
            r = false;
            return
        end
    end
end

r = true;
