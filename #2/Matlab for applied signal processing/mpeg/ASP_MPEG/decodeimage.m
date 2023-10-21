function [im, codedmessage] = decodeimage(codedmessage, imsize, step, varargin)

% DECODEIMAGE Image zigzag/DCT decoder.
%
%   [I, CODEDMESSAGE] = DECODEIMAGE(CODEDMESSAGE, ISIZE, STEP) decodes the char
%   array with the coded DC and AC coeficients of Is' discrete cosine
%   transform, and returns the corresponding image. I is a N-by-M (multiple
%   of 8) uint8 matrix.
%
%   [I, CODEDMESSAGE] = DECODEIMAGE(CODEDMESSAGE, ISIZE, STEP, TYPE) decodes the
%   char array with the coded DC and AC coeficients depending on TYPE:
%
%       TYPE =        
%               'image' - Does the same as above.
%
%               'difference' - CODEDMESSAGE contains the coded DC and AC
%               coeficients of Is' discrete cosine transform, and returns
%               the corresponding differences between two images.
%
%               'motiondifference' - Its similar to 'difference' but uses
%               other huffman codes.
%
%               In the last two cases I is a N-by-M, double matrix that
%               represents the diference between two images.
%
% See ENCODEIMAGE.

if nargin > 3
    type = varargin{1};
else
    type = 'image';
end
    
% Calculate the correspondig zigzag dct size, dcsize and number of DC
% and AC coeficients.
zzsize = imsize .* [1/8 8];
dcsize = [zzsize(1, 1) ceil(zzsize(1, 2) / 64)];
ndc = prod(dcsize);
nac = ndc * 63;

% Extract header from codedmessage and calculate number of ACNZ and ACVV to
% decode.
nbits = ceil(log2(nac));
nzlen = bin2dec(codedmessage(1, 1 : nbits));
vvlen = bin2dec(codedmessage(1, nbits + 1 : 2 * nbits));
codedmessage = codedmessage(1, 2 * nbits + 1 : end);

switch type
    case 'image'
        decoderpath = './decoders/im/';
    case 'difference'
        decoderpath = './decoders/imdiff/';
    case 'motiondifference'
        decoderpath = './decoders/immotiondiff/';
    otherwise
        error('Type not supported!! Check help decodeimage!!');
end

load(strcat(decoderpath, 'dcddecoder.mat'));
[dcdindex codedmessage] = decodemessage(codedmessage, dcddecoder, ndc);
clear dcddecoder;

load(strcat(decoderpath, 'nzdecoder.mat'));
[nzindex codedmessage] = decodemessage(codedmessage, nzdecoder, nzlen);
clear nzdecoder;

load(strcat(decoderpath, 'vvdecoder.mat'));
[vvindex codedmessage] = decodemessage(codedmessage, vvdecoder, vvlen);
clear vvdecoder;

% Calculate DCT - DC components from DC differentiated.
dcd = reshape(dcdindex, dcsize) - 256;
dc = zeros(dcsize + [1 2]);
for i = 2 : size(dc,1)
    for j = 2 : size(dc,2) - 1
        dc(i,j) = ceil(dcd(i-1,j-1) + 1/4 * (dc(i, j-1) + dc(i-1, j-1) + dc(i-1, j) + dc(i-1, j+1)));
    end
end
dc = dc(2 : end, 2 : end - 1);

% Calculate AC coeficients from AC number of zeros (nz) and AC
% values (vv).
ac = irunlength(nzindex - 2, vvindex - 129, 63, nac);
ac = reshape(ac, [nac/zzsize(1, 1) zzsize(1, 1)])' + 128;

% Join DC and AC components into a new DCT.
zzdctq = zeros(zzsize);
mask = true(zzsize);
mask(:, 1 : 64 : end) = false;
zzdctq(~mask) = dc;
zzdctq(mask) = ac(:);

% Calculate image using the DCT calculated above.
if strcmp(type, 'image')
    im = uint8(blkproc(zzdctq, [1 64], @izigzagdctq8, step));
else
    im = round(blkproc(zzdctq, [1 64], @izigzagdctq8diff, step));
end

function im = izigzagdctq8(dctqd, step)
    im = idct2(iquantizedct8(izigzag8(dctqd), step));

function imdiff = izigzagdctq8diff(dctqd, step)
    imdiff = idct2(iquantizedct8diff(izigzag8(dctqd), step));