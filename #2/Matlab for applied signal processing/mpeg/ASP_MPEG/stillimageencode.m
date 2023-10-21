function [Ir,bits,r] = stillimageencode(I,k)

% Performs JPEG coding.
%

if (nargin<2),
    k = 1;
end;

% Build file header: 5 bytes (40 bits) header specifing image size (9 bits
% for each dimension - max 512x512
header(1,  1 : 20) = '0';
header(1,  1 : 10) = dec2bin(size(I,1),10);
header(1, 11 : 20) = dec2bin(size(I,2),10);

% Factor used to operate the diferent Y' and Cb', Cr' picture componentes.
%cf = [16 8 8];

% Steps used in DCT quantization.
step1 = k*jpegsteps;
step2 = k*jpegsteps; step2(1, 1) = step2(1, 1) * 2;

ycbcr = RGBtoYCbCrimage(I,[2 2]);

% Build payload segment.
cmsize = 100000;
codedframesmessage(1, 1 : cmsize) = '0';
codedframesmessage(1, 1 : 20) = header(1, 1 : end);
cmfilled = 20;

[codedframesmessage,cmsize,cmfilled] = channelencode(ycbcr.Y,step1,codedframesmessage,cmsize,cmfilled);
[codedframesmessage,cmsize,cmfilled] = channelencode(ycbcr.Cb,step1,codedframesmessage,cmsize,cmfilled);
[codedframesmessage,cmsize,cmfilled] = channelencode(ycbcr.Cr,step1,codedframesmessage,cmsize,cmfilled);

% PACKBITS! AND CUT OF UNEEDED CODEFRAMESMESSAGE SPACE!!!
%codedframes = packbits(codedframesmessage(1, 1 : cmfilled));
%bits = codedframes;
%codedframes = packbits(codedframesmessage(1, 1 : cmfilled));
bits = codedframesmessage(1,1:cmfilled);

% Reconstructed image
ycbcr.Y  = blkproc(ycbcr.Y, [8 8], @dct_q_iq_idct8, step1);
ycbcr.Cb = blkproc(ycbcr.Cb, [8 8], @dct_q_iq_idct8, step1);
ycbcr.Cr = blkproc(ycbcr.Cr, [8 8], @dct_q_iq_idct8, step1);

Ir = YCbCrtoRGBimage(ycbcr);

% Psnr
r = mypsnr(I,Ir);


function [codedframesmessage,cmsize,cmfilled] = channelencode(I,step1,codedframesmessage,cmsize,cmfilled)


% Encode one picture.
Icodedmessage = encodeimage(applypadding(I,8),step1);
                
% Write coded message in the payload segment.
len = length(Icodedmessage);
while(len > cmsize - cmfilled),
    codedframesmessage(1, end + 1 : end + cmsize) = '0';
    cmsize = 2 * cmsize;
end
codedframesmessage(1, cmfilled + 1 : cmfilled + len) = Icodedmessage(1, :);
cmfilled = cmfilled + len;

                   



function im = applypadding(im, max)
padding = mod(max - mod(size(im), max), max);
im = im([1 : end end - 1 : -1 : end - padding(1, 1)], ...
    [1 : end end - 1 : -1 : end - padding(1, 2)]);

function imr = zigzagdctq8(im, step)
imr = uint8(idct2(iquantizedct8(quantizedct8(dct2(im), step), step)));

