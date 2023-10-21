function RGB709 = YCbCrtoRGB709gc(YCbCr, option)

%YCbCrtoRGB709gc   Conversion of ITU-R BT.601-2 Y'CbCr to gamma corrected ITU-R 
%                  BT.709-2 R'G'B'.
%
%   RGB709 = YCbCrtoRGB709gc(YCbCr) returns convertion from Recommendation
%   ITU-R BT.601-2 Y'CbCr to Recommendation ITU-R BT.709-2 gamma corrected 
%   R'G'B'. 
%  
%   YCbCr must be a uint8 array with exactly three lines, corresponding to
%   Y', Cb, and Cr.
%
%   The R', G', and B' values are returned in an uint8 array of the same size as 
%   YCbCr and are normalized to the interval [0, 255].
%
%   RGB709 = YCbCrtoRGB709gc(YCbCr, 'double') returns a double array with the
%   R', G', and B' values normalized to the interval [0, 1].
%
%   Transformation taken from "Frequently Asked questions about Color", by 
%   Charles Poynton (see http://www.poynton.com/ColorFAQ.html).
%
%   See also RGB709gctoYCbCr.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if ~isa(YCbCr, 'uint8')
    error('YCbCr must be uint8 array');
elseif size(YCbCr, 1) ~= 3
    error('YCbCr must have three lines');
elseif nargin == 2 & ~strcmp(option, 'double')
    error('option must be "double"');
end

s = size(YCbCr);

YCbCr = YCbCr(:, :);

M = [ 65.481 128.553 24.966
     -37.797 -74.203 112
     112 -93.786 -18.214];

if nargin ~= 2
    RGB709 = reshape(uint8(limit(M \ (double(YCbCr) - ...
        [16; 128; 128] * ones(1, size(YCbCr, 2))) * 255, 0, 255)), s);
else
    RGB709 = reshape(limit(M \ (double(YCbCr) - ...
        [16; 128; 128] * ones(1, size(YCbCr, 2))), 0, 1), s);
end
