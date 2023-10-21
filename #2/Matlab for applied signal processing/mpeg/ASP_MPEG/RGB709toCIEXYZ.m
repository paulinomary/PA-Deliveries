function CIEXYZ = RGB709toCIEXYZ(RGB709)

%RGB709toCIEXYZ  Conversion of ITU-R BT.709-2 RGB to CIE XYZ.
%
%   CIEXYZ = RGB709toCIEXYZ(RGB709) returns the conversion of tristimulus values
%   from RGB according to Recommendation ITU-R BT.709-2 to CIE XYZ.  The 
%
%   The BT.709-2 tristimulus values R, G and B in RGB709 must be normalized to
%   the interval [0, 1].  The white point corresponds to R = 1, G = 1, and B = 1.
%   Hence, the tristimulus values are dimensionless.  
%
%   The array RGB709 must have exactly three lines, corresponding to R, G, and B.
%
%   The resulting CIE XYZ tristimulus values are always valid.
%
%   Transformation taken from "Frequently Asked questions about Color", by Charles 
%   Poynton (see http://www.poynton.com/ColorFAQ.html).
%
%   See also CIEXYZtoRGB709.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if ~isa(RGB709, 'double')
    error('RGB709 must be a double array');
elseif size(RGB709, 1) ~= 3
    error('RGB709 must have three lines');
end

s = size(RGB709);

RGB709 = RGB709(:, :);

M = [0.412453 0.357580 0.180423
     0.212671 0.715160 0.072169
     0.019334 0.119193 0.950227];

CIEXYZ = reshape(M * RGB709, s);
