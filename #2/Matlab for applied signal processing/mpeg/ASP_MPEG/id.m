function varargout = id(varargin)

%id   Identity.
%
%   [a, b, c, ...] = id(x, y, z, ...) simply matches up the input and
%   output lists.  It is the same as a = x, b = y, c = z, ...
%
%   Nice idiom for swaping values:
%      [a, b] = id(b, a);
%
%   See also deal.
%
%Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

if nargin ~= nargout
    error('Input and output size must match');
end

varargout = varargin;