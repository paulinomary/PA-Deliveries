function RGB709 = YtoRGB709gc(Y, option)

%YtoRGB709gc   Conversion of ITU-R BT.601-2 Y' to gamma corrected ITU-R 
%                  BT.709-2 R'G'B'.
%
%   RGB709 = YtoRGB709gc(Y) returns convertion from Recommendation
%   ITU-R BT.601-2 Y' to Recommendation ITU-R BT.709-2 gamma corrected 
%   R'G'B'. 
%  
%   'Y' must be a uint8 array with exactly one line, corresponding to Y'.
%
%   The R', G', and B' values are returned in an uint8 array of the same size as 
%   'Y' and are normalized to the interval [0, 255].
%
%   RGB709 = YtoRGB709gc(Y, 'double') returns a double array with the
%   R', G', and B' values normalized to the interval [0, 1].
%
%   Transformation taken from "Frequently Asked questions about Color", by 
%   Charles Poynton (see http://www.inforamp.net/~poynton/).
%
%   See also RGB709gctoY, RGB709gctoYCbCr, YCbCrtoRGB709gc.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if ~isa(Y, 'uint8')
    error('Y must be uint8 array');
elseif size(Y, 1) ~= 1
    error('YCbCr must have one line');
elseif nargin == 2 & ~strcmp(option, 'double')
    error('option must be "double"');
end

s = size(Y);
s(1) = 3;

Y = Y(:, :);

M = [ 65.481 128.553 24.966
     -37.797 -74.203 112
     112 -93.786 -18.214];

M = inv(M);

M = M(:, 1);

if nargin ~= 2
    RGB709 = reshape(uint8(limit(M * (double(Y) - ...
        16 * ones(1, size(Y, 2))) * 255, 0, 255)), s);
else
    RGB709 = reshape(limit(M * (double(Y) - ...
        16 * ones(1, size(Y, 2))), 0, 1), s);
end
