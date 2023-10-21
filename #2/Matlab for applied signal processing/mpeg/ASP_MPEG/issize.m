function r = issize(s)

% issize - True for size row vectors.
%
%   r = issize(s)  Returns true if 's' is a size row vector, that is, a double 
%   row vector of non-negative values.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

r = isrow(s) & isa(s, 'double') & all(0 <= s);