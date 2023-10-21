function RGB = YCbCrtoRGBimage(YCbCr, varargin)

%YCbCrtoRGBimage   Conversion from Y'CbCr to R'G'B' Rec. 709.
%
%   RGB = YCbCrtoRGBimage(YCbCr, ...) returns an image array where each image
%   is the conversion of one image in struct array 'YCbCr' to R'G'B'.  Each 
%   element of the struct array corresponds to one R'G'B' image in the returned
%   array and consists of four fields:
%      Y - the Y' component.
%      Cb - the Cb component.
%      Cr - the Cr component.
%      ss - the same value as 'ss' (sub-sampling factor).
%   The chrominance components in YCbCr may be sub-sampled.  The sub-sampling
%   factors are obtained by comparing the sizes of the chrominances to
%   the size of the luminance (which is not sub-sampled).  However, the field
%   'ss' must be present in 'YCbCr'.
%
%   Extra arguments are passed to imresize (e.g., 'bilin').
%
%   Example:
%
%    Read an image:
%
%      balloonRGB = readimage('.../balloon.bmp');
%
%    Convert from R'G'B' to Y'CbCr, subsampling chrominances by 8 in each 
%    direction:
%
%      balloonYCbCr = RGBtoYCbCrimage(balloonRGB, [8 8], 'bilin');
%
%    Show luminance component (stretched to fill [0, 255] interval):
%
%      showimage(strecthimage(balloonYCbCr.Y));
%
%    Convert back to R'G'B':
%
%      balloonrec = YCbCrtoRGBimage(balloonYCbCr, 'bilin');
%
%    Compare to orginal:
%
%      showimage(balloon);
%      showimage(balloonrec);
%
%   See imresize, RGBtoYCbCrimage, RGBtoYimage, YtoRGBimage.
%
%Copyright (c) 2002 - Manuel Menezes de Sequeira (ISCTE).


% Verificar
% YCbCr - ????

% 1. Upsample

cs = size(YCbCr);
i = ones(size(cs));
for j = 1 : prod(cs)
    ci = num2cell(i);
    ss = YCbCr(ci{:}).ss;
    s = size(YCbCr(ci{:}).Y);
    if any(ss ~= [1 1])
        YCbCr(ci{:}).Cb = uint8(128 + ...
            imresize(double(YCbCr(ci{:}).Cb) - 128, s, varargin{:}));
        YCbCr(ci{:}).Cr = uint8(128 + ...
            imresize(double(YCbCr(ci{:}).Cr) - 128, s, varargin{:}));
        i = next(i - 1, cs) + 1;
    end
end

% 2. Concatenate into one, big, array:

Y = cat(4, YCbCr.Y);
Cb = cat(4, YCbCr.Cb);
Cr = cat(4, YCbCr.Cr);

YCbCr = cat(3, Y, Cb, Cr);

YCbCr = reshape(YCbCr, [s 3 cs]);

RGB = ipermute(YCbCrtoRGB709gc(permute(YCbCr, [3 1 2 4 : ndims(YCbCr)])),...
    [3 1 2 4 : ndims(YCbCr)]);

return

% 2. Non-cell,permute,


is = size(RGB);

im = RGB709gctoYCbCr(permute(RGB, [3 1 2 4 : ndims(RGB)]));

iss = size(im);
iss(1) = 1;
Y = num2cell(ipermute(reshape(im(1, :), iss), [3 1 2 4 : ndims(RGB)]), [1 2]);
Cb = num2cell(ipermute(reshape(im(2, :), iss), [3 1 2 4 : ndims(RGB)]), [1 2]);
Cr = num2cell(ipermute(reshape(im(3, :), iss), [3 1 2 4 : ndims(RGB)]), [1 2]);


if length(is) < 4
    is(4) = 1;
end

YCbCr(is(4 : end)).Y = 0;
YCbCr(is(4 : end)).Cb = 0;
YCbCr(is(4 : end)).Cr = 0;

[YCbCr(:).Y] = id(Y{:});
[YCbCr(:).Cb] = deal(Cb{:});
[YCbCr(:).Cr] = deal(Cr{:});

YCbCr.ss = ss;