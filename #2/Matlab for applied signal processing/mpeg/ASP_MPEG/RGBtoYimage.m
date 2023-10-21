function Y = RGBtoYimage(RGB)

%RGBtoYimage   Conversion from R'G'B' Rec. 709 to Y'.
%
%   Y = RGBtoYimage(RGB) returns an array of images with the same dimensions
%   as 'RGB' where each image is the conversion of one image in 'RGB'
%   to Y'.  The returned array has the same size as 'RGB' except for the third
%   dimension, whose size will be 1 instead of 3.
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
%   Notice access using dot ('.') operator.
%
%   See imresize, RGBtoYCbCrimage, YCbCrtoRGBimage, YtoRGBimage.
%
%Copyright (c) 2002 - Manuel Menezes de Sequeira (ISCTE).

% Atenç~ao!  A imagem devolvida n~ao usa toda a gama dinamica!  Ver Rec. 601.
% Para mostrar no ecra usar showimage(stretchimage(xxx)).
% Verificar
% ss length 2, ss >= 1, row, float
% RGB - uint8 or double

Y = ipermute(RGB709gctoY(permute(RGB, [3 1 2 4 : ndims(RGB)])),...
    [3 1 2 4 : ndims(RGB)]);
