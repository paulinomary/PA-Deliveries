function [maxind,minind] = findpeaks(x)

% maxind = findpeaks(x) returns the indexes of the local maxima in array x.
% [maxind,minind] = findpeaks(x) also returns the indexes of the minima

maxind = 1+find(x(2:end-1) > x(3:end) & x(2:end-1) >= x(1:end-2));
if(nargout > 1 | nargin > 1)
   minind = 1+find(x(2:end-1) < x(3:end) & x(2:end-1) <= x(1:end-2));
end
