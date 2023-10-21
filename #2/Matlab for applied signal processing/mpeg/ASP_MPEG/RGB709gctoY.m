function Y = RGB709gctoY(RGB709)

%RGB709gctoY   Conversion of gamma corrected ITU-R BT.709-2 R'G'B' to ITU-R
%                  BT.601-2 Y'.
%
%   Y = RGB709gctoY(RGB) returns the conversion from ITU-R Recommendation
%   BT.709-2 gamma corrected R'G'B' to ITU-R Recommendation BT.601-2 Y'.  
%
%   The R', G' and B' values in RGB709 must be normalized to the interval [0, 1]
%   for double arrays and to the interval [0, 255] for uint8 arrays. 
%
%   RGB709 must be an array with exactly three lines, corresponding to R', G', 
%   and B'.
%
%   Notice that the Y' values are uint8 and do not ocupy the whole [0, 255] 
%   range.
%
%   Transformation taken from "Frequently Asked questions about Color", by 
%   Charles Poynton (see http://www.inforamp.net/~poynton/).
%
%   See also YtoRGB709gc, RGB709gctoYCbCr, YCbCrtoRGB709gc.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if ~isa(RGB709, 'uint8') & ~isa(RGB709, 'double')
    error('RGB709 must be an uint8 or double array');
elseif size(RGB709, 1) ~= 3
    error('RGB709 must have three lines');
end

s = size(RGB709);
s(1) = 1;

RGB709 = RGB709(:, :);

M = [ 65.481 128.553 24.966];
if isa(RGB709, 'uint8')
    Y = reshape(uint8(limit(16 * ones(1, size(RGB709, 2)) +...
        M / 255 * double(RGB709), 0, 255) + 0.5), s);
else
    Y = reshape(uint8(limit(16 * ones(1, size(RGB709, 2)) +...
        M * RGB709, 0, 255) + 0.5), s);
end