function packed = slowpackbits(bits)

%slowpackbits   Packs a binary string into a uint8 row vector.
%
%   packed = packbits(bits) returns a uint8 row vector containing the same 
%   information as the binary string 'bits'.  Zero bits are stufed until the 
%   the length is a multiple of 8.
%
%   It is slower than 'packbits', though it deals better with large amounts 
%   of information.
%
%   See packbits.
% 
%Copyright (C) 2002 - Manuel Menezes de Sequeira (ISCTE).

if(nargin ~= 1)
    error('must have 1 argument: bits to pack');
end

if(~isbinarychar(bits) | ~isrow(bits))
    error('bits must be a single binary string (row vector)');
end
        
% Stuff with zeros until a multiple of 8 bits is obtained:
extra = 8 - mod(length(bits), 8);
if(extra ~= 8)
    bits(end + 1 : end + extra) = '0';
end

bits = reshape(bits, [8 length(bits) / 8]);

packed = uint8(logical(bits(8, :)) + ...
    logical(bits(7, :)) * 2 + ...
    logical(bits(6, :)) * 4 + ...
    logical(bits(5, :)) * 8 + ...
    logical(bits(4, :)) * 16 + ...
    logical(bits(3, :)) * 32 + ...
    logical(bits(2, :)) * 64 + ...
    logical(bits(1, :)) * 128);
