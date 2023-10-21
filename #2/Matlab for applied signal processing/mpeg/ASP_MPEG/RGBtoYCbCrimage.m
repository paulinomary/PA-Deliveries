function YCbCr = RGBtoYCbCrimage(RGB, ss, varargin)

%RGBtoYCbCrimage   Conversion from R'G'B' Rec. 709 to Y'CbCr.
%
%   YCbCr = RGBtoYCbCrimage(RGB, ss, ...) returns a struct array where each
%   element is the conversion of one image in 'RGB' to Y'CbCr.  Each element of
%   the struct array corresponds to one R'G'B' image in 'RGB' and consists of
%   three fields:
%      Y - the Y' component.
%      Cb - the Cb component.
%      Cr - the Cr component.
%      ss - the same value as 'ss' (sub-sampling factor).
%   The chrominance components are sub-sampled acording to 'ss' (e.g. [2 4]
%   means sub-sample verticaly by 2 and horizontaly by 4).
%
%   The Y', Cb, and Cr matrices in the structure are uint8 and do not ocupy the
%   whole [0, 255] range!
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
%      showimage(stretchimage(balloonYCbCr.Y));
%
%   Notice access using dot ('.') operator.
%
%   See imresize, YCbCrtoRGBimage, RGBtoYimage, YtoRGBimage.
%
%Copyright (c) 2002 - Manuel Menezes de Sequeira (ISCTE).

% Obsolete?  See rgb2ycbcr!  Mention rgb2gray as a similar transformation
% which fills the whole [0, 255] range.

% Verificar
% ss length 2, ss >= 1, row, float
% RGB - uint8 or double

if nargin < 2
    ss = [1 1];
end

is = size(RGB);

im = RGB709gctoYCbCr(permute(RGB, [3 1 2 4 : ndims(RGB)]));

iss = size(im);
iss(1) = 1;
Y = num2cell(ipermute(reshape(im(1, :), iss), [3 1 2 4 : ndims(RGB)]), [1 2]);
Cb = num2cell(ipermute(reshape(im(2, :), iss), [3 1 2 4 : ndims(RGB)]), [1 2]);
Cr = num2cell(ipermute(reshape(im(3, :), iss), [3 1 2 4 : ndims(RGB)]), [1 2]);

if any(ss ~= [1 1])
    cs = size(Cb);
    i = ones(size(cs));
    for j = 1 : prod(cs)
        ci = num2cell(i);
        Cb{ci{:}} = uint8(128 + ...
            imresize(double(Cb{ci{:}}) - 128, is(1 : 2) ./ ss, varargin{:}));
        Cr{ci{:}} = uint8(128 + ...
            imresize(double(Cr{ci{:}}) - 128, is(1 : 2) ./ ss, varargin{:}));
        i = next(i - 1, cs) + 1;
    end
end

if length(is) < 4
    is(4) = 1;
end

if length(is) < 5
    is(5) = 1;
end

is = num2cell(is(4 : end));

YCbCr(is{:}).Y = 0;
YCbCr(is{:}).Cb = 0;
YCbCr(is{:}).Cr = 0;
YCbCr(is{:}).ss = 0;

[YCbCr.Y] = id(Y{:});
[YCbCr.Cb] = id(Cb{:});
[YCbCr.Cr] = id(Cr{:});

[YCbCr.ss] = deal(ss);