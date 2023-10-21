function packed = packbits(bits)

%packbits   Packs a logical row vector into a uint8 row vector.
%
%   packed = packbits(bits) returns a uint8 row vector containing the same 
%   information as the logical row vector 'bits'.  Zero bits are stufed
%   until the length is a multiple of 8.
% 
%   See slowpackbits.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if(nargin ~= 1)
    error('must have 1 argument: bits to pack');
end

if(~islogical(bits) | ~isrow(bits))
    error('bits must be a logical row vector');
end
        
% Stuff with zeros until a multiple of 8 bits is obtained:
extra = 8 - mod(length(bits), 8);
if(extra ~= 8)
    bits(end + 1 : end + extra) = false;
end

bits8 = double(reshape(bits ~= 0, [8 length(bits) / 8]));

packed = uint8(2 .^ [7 : -1 : 0] * bits8);