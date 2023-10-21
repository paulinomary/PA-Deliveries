function xyz = XYZtoxyz(XYZ)

%XYZtoxyz   Conversion from tristimulus values to chromaticity coordinates.
%
%   xyz = XYZtoxyz(XYZ) returns the conversion of tristimulus values XYZ to 
%   chromaticity coordinates.
%
%   XYZ must be an array an array with exactly three lines, corresponding to X,
%   Y, and Z.
%
%   Returns an array of the same size as XYZ with the x, y, and z chromaticity
%   coordinates.
%   
%   The conversion goes as follows:
% 
%       x = X / (X + Y + Z);
%       y = Y / (X + Y + Z);
%       z = Z / (X + Y + Z);
% 
%   This transformation corresponds to the central projection of the tristimulus
%   values onto a plane which intersects the three axes x, y and at 1.
%   The values of x, y, and z are thus dependent.  The first two are usualy taken
%   to draw chromaticity diagrams.  The third may be calculated from them by:
%
%       z = 1 - x - y;
%
%   Since the tristimulus values X, Y, and Z are divided by the their sum, their
%   absolute value is irrelevant.  Only the relative values matter.  Hence, 
%   brightness is absent from the chromaticity coordinates.
%
%   The tristimulus values are usualy taken to lie in the first octant of the 3D
%   space.  This means that they correspond to a color which is physicaly 
%   realizable with their corresponding primaries (it must be remembered that
%   the primaries themselves may be imaginary: that is the case with CIE XYZ).
%   Hence, the resulting chromaticity coordinates usualy have non-negative 
%   values.  
%
%   The XYZ tristimulus does not have to correspond to CIE XYZ.  It must be
%   tristimulus values, but of an arbitrary set of three linearly independent 
%   primaries (also, no linear combination of the primaries can belong to the
%   null space of the cone sensitivities).
%
%   Usualy the tristimulus values are given in cd/m^2 (candelas per square 
%   meter).  The chromaticity corrdinates are dimensionless.
%
%   The "white color" point is taken to be the result of tranforming [1; 1; 1],
%   i.e., the tristimulus values X = 1, Y = 1, and Z = 1.  Hence, it lies on 
%   [1/3, 1/3; 1/3].
%
%   See "Frequently Asked questions about Color", by Charles Poynton (see
%   http://www.inforamp.net/~poynton/).
%
%   See also xyYtoXYZ.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

xyz(1,:) = XYZ(1,:)./(XYZ(1,:)+XYZ(2,:)+XYZ(3,:));
xyz(2,:) = XYZ(2,:)./(XYZ(1,:)+XYZ(2,:)+XYZ(3,:));
xyz(3,:) = XYZ(3,:)./(XYZ(1,:)+XYZ(2,:)+XYZ(3,:));
