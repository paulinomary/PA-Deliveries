function r = logical(s)

% logical - converts a string into a logical matrix.  
% 
%   r = logical(s) returns a double matrix the same size of s with logical value
%   true where 's' is '1' and false where it is '0'.  It is an error if 's' is
%   not a binary char array (binary string).
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).


if(nargin ~= 1)
    error('logical requires one argument');
elseif(~isbinarychar(s))
    error('argument must be binary char array');
end

r = logical(s - '0');
% Should be 'r = s ~= '0';', but ~= behaves badly with empty arguments. 
