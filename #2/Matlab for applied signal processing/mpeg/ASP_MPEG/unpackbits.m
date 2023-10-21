function bits = unpackbits(packed, format)

%unpackbits   Unpacks a uint8 row vector into a binary string.
%
%
%   bits = unpackbits(packed [, format]) returns a
%       a) binary string - if 'format' is 'char';
%       b) logical uint8 row vector - if 'format' is 'uint8'; and
%       c) logical double row vector - if 'format' is 'double'
%   containing the same information as the uint8 row vector 'packed'.  By default
%   format is 'char'.
%
%   See also packbits and slowpackbits.
% 
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if(nargin ~= 1 & nargin ~= 2)
    error('must have 1 or 2 arguments: uint8 packed bits to unpack [and format]');
end

if(~isa(packed, 'uint8') | ~isrow(packed))
    error('packed must be a uint8 row vector');
end

if nargin < 2
    format = 'char';
end
        
table = dec2bin(0 : 255)';

switch lower(format)
case 'char'
case 'uint8'
    table = logical(uint8(table ~= '0'));
case 'double'
    table = logical(table ~= '0');
end

bits = table(:, double(packed) + 1);
bits = bits(:)';