function [codemessage,imQ] = encodeimage(im, stepQ, varargin)

% ENCODEIMAGE Image zigzag/DCT encoder.
%
%   C = ENCODEIMAGE(I, stepQ) returns a char row vector with the coded DC
%   and AC coeficients of Is' discrete cosine transform. I must be a
%   N-by-M uint8 matrix with size multiple of 8. stepQ is a 8-by-8 matrix
%   that contains the quantization steps.
%
%   C = ENCODEIMAGE(I, stepQ, TYPE) type may be specified:
%
%       TYPE =        
%               'image' - Does the same as above.
%
%               'difference' - Does the same as above but since I is the
%               difference between to images, uses another codification
%               method.
%
%               'motiondifference' - Its similar to 'difference' but in
%               this case the difference is between an image estimated with
%               motion vectors and the original one.
%
% See DECODEIMAGE

if nargin > 2        
    type = varargin{1};
else
    type = 'image';
end

imsize = size(im);

% Get DCT components from 8x8 blocks of the image, in the form of 1x64
% zigzag blocks.
switch type
    case 'image'
        if isa(im, 'uint8')
            zzdctq = double(blkproc(im, [8 8], @zigzagdctq8, stepQ));
            load('./codes/imcodes.mat')
        else
            error('Image must be a uint8 matrix!');
        end
    case 'difference'
        if isa(im, 'double')
            zzdctq = double(blkproc(im, [8 8], @zigzagdctq8diff, stepQ));
            load('./codes/imdiffcodes.mat')
        else
            error('Image must be a double matrix!');
        end
    case 'motiondifference'
        if isa(im, 'double')
            zzdctq = double(blkproc(im, [8 8], @zigzagdctq8diff, stepQ));
            load('./codes/immotiondiffcodes.mat')
        else
            error('Image must be a double matrix!');
        end
    otherwise
        error('Type is not supported!! Check help encodeimage!!')
end

imQ = blkproc(zzdctq,[1 64],@izigzag8);

[zznrows zzncolumns] = size(zzdctq);

% Separate DC from AC components of the DCT.
mask = true([zznrows zzncolumns]);
mask(:, 1 : 64 : end) = false;
dc = zzdctq(~mask);
ac = reshape(zzdctq(mask), ...
    [zznrows, zzncolumns - zzncolumns / 64])' - 128;

% The DC component level is obtain by the mean difference between is own
% value and the above and left coeficients.
h = [-0.25 -0.25 -0.25; -0.25 1 0; 0 0 0];
dc = reshape(dc, imsize / 8);
dcdindex = floor(filter2(h, dc)) + 256;

% Get NZ and VV sequence of AC component.
[nz vv] = runlength(ac(:)', 63);

% Encode DC, ACNZ and ACVV components with the huffman code.
codeddcd = encodemessage(dcdindex(:)', dcdcode);
codednz = encodemessage(nz(:)' + 2, nzcode);
codedvv = encodemessage(vv(:)' + 129, vvcode);

% Header with codedacnz and codedacvv length
nbits = ceil(log2(63 * prod(imsize) / 64));
header = strcat(dec2bin(length(nz), nbits), dec2bin(length(vv), nbits));

codemessage = strcat(header, codeddcd, codednz, codedvv);

function dctqd = zigzagdctq8(im, stepQ)
dctqd = zigzag8(quantizedct8(dct2(im), stepQ));

function dctqd = zigzagdctq8diff(imdiff, stepQ)
dctqd = zigzag8(quantizedct8diff(dct2(imdiff), stepQ));
