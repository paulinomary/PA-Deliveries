function RGB = YtoRGBimage(Y)

%YtoRGBimage   Conversion from Y' to R'G'B' Rec. 709.
%
%   RGB = YtoRGBimage(Y) returns an array of images with the same dimensions
%   as 'Y' where each image is the conversion of one image in 'Y'
%   to R'G'B'.  The returned array has the same size as 'Y' except for the third
%   dimension, whose size will be 3 instead of 1.  The resulting image will 
%   display as grayscale, even though it is stored as R'G'B'.
%
%   Example:
%
%    Read an image:
%
%      balloonRGB = readimage('.../balloon.bmp');
%
%    Convert from R'G'B' to Y':
%
%      balloonY = RGBtoYimage(balloonRGB);
%
%    Show obtained luminance component (stretched to fill [0, 255] interval):
%
%      showimage(strecthimage(balloonY));
%
%    Convert to R'G'B' (gray) image:
%
%      balloonrec = YtoRGBimage(balloonY);
%
%    Show result (no stretching necessary):
%
%      showimage(balloonrec);
%
%   Notice access using dot ('.') operator.
%
%   See imresize, RGBtoYCbCrimage, YCbCrtoRGBimage, RGBtoYimage.
%
%Copyright (c) 2002 - Alguem (ISCTE).

RGB = ipermute(YtoRGB709gc(permute(Y, [3 1 2 4 : ndims(Y)])),...
    [3 1 2 4 : ndims(Y)]);