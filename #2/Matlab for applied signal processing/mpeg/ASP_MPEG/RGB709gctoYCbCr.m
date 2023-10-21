function YCbCr = RGB709gctoYCbCr(RGB709)

%RGB709gctoYCbCr   Conversion of gamma corrected ITU-R BT.709-2 R'G'B' to ITU-R
%                  BT.601-2 Y'CbCr.
%
%   YCbCr = RGB709gctoYCbCr(RGB) returns the conversion from ITU-R Recommendation
%   BT.709-2 gamma corrected R'G'B' to ITU-R Recommendation BT.601-2 Y'CbCr.  
%
%   The R', G' and B' values in RGB709 must be normalized to the interval [0, 1]
%   for double arrays and to the interval [0, 255] for uint8 arrays. 
%
%   RGB709 must be an array with exactly three lines, corresponding to R', G', 
%   and B'.
%
%   Notice that the Y'CbCr values are uint8 and do not ocupy the whole [0, 255] 
%   range.
%
%   Transformation taken from "Frequently Asked questions about Color", by 
%   Charles Poynton (see http://www.poynton.com/ColorFAQ.html).
%
%   See also YCbCrtoRGB709gc.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if ~isa(RGB709, 'uint8') & ~isa(RGB709, 'double')
    error('RGB709 must be an uint8 or double array');
elseif size(RGB709, 1) ~= 3
    error('RGB709 must have three lines');
end

s = size(RGB709);

RGB709 = RGB709(:, :);

M = [ 65.481 128.553 24.966
    -37.797 -74.203 112
    112 -93.786 -18.214];

if isa(RGB709, 'uint8')
    YCbCr = reshape(uint8(limit([16; 128; 128] * ones(1, size(RGB709, 2)) +...
        M / 255 * double(RGB709), 0, 255) + 0.5), s);
else
    YCbCr = reshape(uint8(limit([16; 128; 128] * ones(1, size(RGB709, 2)) +...
        M * RGB709, 0, 255) + 0.5), s);
end