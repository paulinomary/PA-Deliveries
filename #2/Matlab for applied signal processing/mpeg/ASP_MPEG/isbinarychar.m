function r = isbinarychar(s)

% isbinarychar - returns true if the argument is a char array containing only
% '0' and '1' characters, i.e., if it is a binary string (matrix).
%
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

if(~ischar(s))
    r = false;
    return
end

% This test should not be necessary, but == behaves badly when on argument is 
% empty...

if(isempty(s))
    r = true;
    return
end

r = s == '0' | s == '1';

r = all(r(:));