function r = limit(m, mini, maxi)

% limit - Returns an array which is equal to its first argument, though with 
% values limited to an interval given by the second an third arguments.
%
% Usage:
%    r = limit(v, min, max)
% where
%    v - array whose values are to be limited.
%    min, max - scalars containing the minimum and maximum allowed values.
%
% Returns the limited version of the input array v.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

%%%%%r = (m > max) .* max + (m < min) .* min + (min <= m & m <= max) .* m;

if fullany(mini >= maxi)
    error('min must be smaller or equal to max');
else
    r = min(max(m, mini), maxi);
end