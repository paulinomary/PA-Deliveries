function [Ir,bits,r] = stillimagegrayencode(I,k)

% JPEG_GRAY_ENCODE Performs JPEG codification.
%

if (nargin<2),
    k = 1;
end;

% Build file header: 5 bytes (40 bits) header specifing image size (9 bits
% for each dimension - max 512x512
header(1,  1 : 20) = '0';
header(1,  1 : 10) = dec2bin(size(I,1),10);
header(1, 11 : 20) = dec2bin(size(I,2),10);

% Build payload segment.
cmsize = 100000;
codedframesmessage(1, 1 : cmsize) = '0';
codedframesmessage(1, 1 : 20) = header(1, 1 : end);
cmfilled = 20;

% Factor used to operate the diferent Y' and Cb', Cr' picture componentes.
%cf = [16 8 8];

% Steps used in DCT quantization.
step1 = k*jpegsteps;
step2 = k*jpegsteps; step2(1, 1) = step2(1, 1) * 2;

% Low-Pass filters
%h = [1 2 1; 2 4 2; 1 2 1] / 16;
%h2 = [0 1 0; 1 4 1; 0 1 0] / 8;

% Encode I picture.
Icodedmessage = encodeimage(applypadding(I,8),step1);
                
% Write coded message in the payload segment.
len = length(Icodedmessage);
while(len > cmsize - cmfilled),
    codedframesmessage(1, end + 1 : end + cmsize) = '0';
    cmsize = 2 * cmsize;
end
codedframesmessage(1, cmfilled + 1 : cmfilled + len) = Icodedmessage(1, :);
cmfilled = cmfilled + len;

                   
% PACKBITS! AND CUT OF UNEEDED CODEFRAMESMESSAGE SPACE!!!
%codedframes = packbits(codedframesmessage(1, 1 : cmfilled));
%bits = codedframes;
bits = codedframesmessage(1,1:cmfilled);

% Reconstructed image
Ir = blkproc(I, [8 8], @dct_q_iq_idct8, step1);

% Psnr
r = psnr(I,Ir);


function im = applypadding(im, max)
padding = mod(max - mod(size(im), max), max);
im = im([1 : end end - 1 : -1 : end - padding(1, 1)], ...
    [1 : end end - 1 : -1 : end - padding(1, 2)]);

function imr = zigzagdctq8(im, step)
imr = uint8(idct2(iquantizedct8(quantizedct8(dct2(im), step), step)));

