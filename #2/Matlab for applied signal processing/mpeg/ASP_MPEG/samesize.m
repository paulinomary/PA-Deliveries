function r = samesize(A, B, ignore)

%samesize   Test equality of size between two variables.
%
%   r = samesize(A, B) returns true if A and B have the same size.
%
%   Two variables have the same size if they have the same size along each 
%   dimension.
%
%   r = samesize(A, B, ignore) returns true if A and B have the same size 
%   ignoring diferences along dimensions given in ignore.  Parameter ignore
%   must be a row vector.
%
%Copyright (C) 2002 - Manuel Menezes de Sequeira (ISCTE).

if nargin < 3
    ignore = zeros(1, 0);
end

if ~isrow(ignore)
    error('ignore must be a row vector');
elseif ~myisinteger(ignore) | any(ignore < 1)
    error('ignore must have only positive integers');
end

sa = size(A);
sb = size(B);

n = max(max(ndims(A), ndims(B)), max(ignore));

if n > ndims(A)
    sa = [sa ones(1, n - ndims(A))];
end

if n > ndims(B)
    sb = [sb ones(1, n - ndims(B))];
end

r = sa == sb;

r(ignore) = [];

r = all(r);