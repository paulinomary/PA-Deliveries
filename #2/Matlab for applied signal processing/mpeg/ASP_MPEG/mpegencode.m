function [codedframes,tnbits,ttype,tpsnr] = mpegencode(codedseqpath, seqpath, first, last, varargin)

% MPEGENCODE Performs MPEG2 codification.
%
%   C = MPEGENCODE(CP, SP, FF, LF)
%       CP - Coded sequency path (output)
%       SP - Sequence path (input)
%       FF - First Frame
%       LF - Last Frame
%
% PP, PI)
%   C = MPEGENCODE(CP, SP, FF, LF, PP, PI) period of P and I frames may be set.
%       PP - Predictive Frame period (B frames in between)
%       PI - Intra Coding Frames period (P frames in between)
%       
% See MPEGDECODE

verbose = 1;
cmfilled_old = 0;

% Set the period of P and I pictures.
pP = 5;
pI = 4;
Qk = 1.0;
if nargin > 4
    pP = varargin{1};
    if nargin > 5
        pI = varargin{2};
        if nargin > 6
            Qk = varargin{3},
        end
    end
end

% Set number of frames.
nframes = last - first + 1;

% Set deviation from 1 to first.
d = first - 1;

% Get sequence size.
seqsize = sequencesize(seqpath);

% Build file header: 5 bytes (40 bits) header specifing image size (9 bits
% for each dimension - max 512x512, number of frames encoded (9 bits - max
% 512 frames), period of P and I pictures (7 bits for P, and 6 bits for
% I bits - max 128 and 64).
header(1, 1 : 40) = '0';
header(1, 1 : 9) = dec2bin(seqsize(1, 1), 9);
header(1, 10 : 18) = dec2bin(seqsize(1, 2), 9);
header(1, 19 : 27) = dec2bin(nframes, 9);
header(1, 28 : 34) = dec2bin(pP, 7);
header(1, 35 : 40) = dec2bin(pI, 6);

% Build payload segment.
cmsize = 100000;
codedframesmessage(1, 1 : cmsize) = '0';
codedframesmessage(1, 1 : 40) = header(1, 1 : end);
cmfilled = 40;

% Set stream codification order
[pictype index] = shufflestream(pP, pI, nframes);

% Factor used to operate the diferent Y' and Cb', Cr' picture componentes.
cf = [16 8 8];

% Sub-sampling factor.
ss = [2 2];

% Load huffmancode used to encode predictive mask.
load('./codes/pcode.mat');

% Steps used in DCT quantization.
step1 = jpegsteps.*Qk;
%step2 = jpegsteps.*Qk; step2(1, 1) = step2(1, 1) * 2;
step2 = 16*ones(size(step1)).*Qk*2;

% Low-Pass filters
h = [1 2 1; 2 4 2; 1 2 1] / 16;
h2 = [0 1 0; 1 4 1; 0 1 0] / 8;

% Y' Cb' Cr' structure.
YCbCrstruct = struct('Y', [], 'Cb', [], 'Cr', [], 'ss', ss);
fnames = fieldnames(YCbCrstruct);

% Backward and foward reference image.
backwardimr = YCbCrstruct;
forwardimr = YCbCrstruct;

for i = 1 : nframes
    
    switch pictype(index(i))
        
        case 'I'
            
            % Update backward image reference.
            backwardimr = forwardimr;
            
            % Get picture
            imI = RGBtoYCbCrimage(...
                readsequence(seqpath, index(i) + d, index(i) + d), ...
                ss, 'bilinear');

            for j = 1 : 3
                
                % Encode I picture.
                Icodedmessage = encodeimage(applypadding( ...
                    getfield(imI, fnames{j}), cf(1, j)), step1);
                
                % Write coded message in the payload segment.
                len = length(Icodedmessage);
                while(len > cmsize - cmfilled)
                    codedframesmessage(1, end + 1 : end + cmsize) = '0';
                    cmsize = 2 * cmsize;
                end
                codedframesmessage(1, cmfilled + 1 : cmfilled + len) = ...
                    Icodedmessage(1, :);
                cmfilled = cmfilled + len;

                % Update forward image reference with decoder results.
                forwardimr = setfield(forwardimr, fnames{j}, ...
                    decoderesult(applypadding( ...
                    getfield(imI, fnames{j}), cf(1 ,j)), step1));
            end

            tnbits(index(i)+d) = cmfilled - cmfilled_old;
            ttype(index(i)+d) = 'I';
            tpsnr.y(index(i)+d) = mypsnr(forwardimr.Y,imI.Y);
            tpsnr.cb(index(i)+d) = mypsnr(forwardimr.Cb,imI.Cb);
            tpsnr.cr(index(i)+d) = mypsnr(forwardimr.Cr,imI.Cr);
            cmfilled_old = cmfilled;            
            if (verbose == 1),
                disp(['Coding ' ttype(index(i)+d) ' frame (' num2str(index(i)+d),'):' ...
                    ' bits=' num2str(tnbits(index(i)+d)) ...
                    ', psnr=[' num2str(tpsnr.y(index(i)+d)) ' ' num2str(tpsnr.cb(index(i)+d)) ' ' num2str(tpsnr.cr(index(i)+d)) ']']);
            end;
            
        case 'P'

            % Update backward image reference.
            backwardimr = forwardimr;

            % Get picture
            imP = RGBtoYCbCrimage( ...
                readsequence(seqpath, index(i) + d, index(i) + d), ...
                ss, 'bilinear');
            
            [mvv mhv errors] = estimatemotion(applypadding(imP.Y, 16), ...
                backwardimr.Y, [16 16], [15 15], 'fastfullsearch');

            % Apply motion estimated to the reference image with and
            imPestimated = applymotion(backwardimr.Y, mvv, mhv, [16 16]);
                           
            % Differences between imP estimated and the original imP.
            imPdiff = double(applypadding(imP.Y, 16)) - ...
                double(imPestimated);
            
            % Encode diferences and motion vectores.
            diffcm = encodeimage(imPdiff, step2, 'motiondifference');
            mvcm = encodemv(mvv, mhv);
                        
            Pcodedmessage = strcat(mvcm, diffcm);
            
            % Write coded message in the payload segment.
            len = length(Pcodedmessage);
            while(len > cmsize - cmfilled)
                codedframesmessage(1, end + 1 : end + cmsize) = '0';
                cmsize = 2 * cmsize;
            end
            codedframesmessage(1, cmfilled + 1 : cmfilled + len) = ...
                Pcodedmessage(1, :);
            cmfilled = cmfilled + len;
            
            % Apply decoder quantization errors to differences.
            imPdiff = decoderesult(imPdiff, step2, 'difference');
            
            % Sum the differences to the image estimated with motion
            % vectors.
            forwardimr.Y = double(imPestimated) + imPdiff;
            
            % Compute and encode Cb' and Cr' components differences.
            for j = 2 : 3

                imPestimated = applymotion( ...
                    getfield(backwardimr, fnames{j}), mvv/2, mhv/2, [8 8]);
                        
                imPdiff = ...
                    double(applypadding(getfield(imP, fnames{j}), 8)) - ...
                    double(imPestimated);
            
                diffcm = encodeimage(imPdiff, step2, 'motiondifference');
                len = length(diffcm);
                while(len > cmsize - cmfilled)
                    codedframesmessage(1, end + 1 : end + cmsize) = '0';
                    cmsize = 2 * cmsize;
                end
                codedframesmessage(1, cmfilled + 1 : cmfilled + len) = ...
                    diffcm(1, :);
                cmfilled = cmfilled + len;

                imPdiff = decoderesult(imPdiff, step2, 'difference');
                forwardimr = setfield(forwardimr, fnames{j}, ...
                    double(imPestimated) + imPdiff);
            end
            
            tnbits(index(i)+d) = cmfilled - cmfilled_old;
            ttype(index(i)+d) = 'P';
            tpsnr.y(index(i)+d) = mypsnr(uint8(forwardimr.Y),imP.Y);
            tpsnr.cb(index(i)+d) = mypsnr(uint8(forwardimr.Cb),imP.Cb);
            tpsnr.cr(index(i)+d) = mypsnr(uint8(forwardimr.Cr),imP.Cr);
            cmfilled_old = cmfilled;            
            if (verbose == 1),
                disp(['Coding ' ttype(index(i)+d) ' frame (' num2str(index(i)+d),'):' ...
                    ' bits=' num2str(tnbits(index(i)+d)) ...
                    ', psnr=[' num2str(tpsnr.y(index(i)+d)) ' ' num2str(tpsnr.cb(index(i)+d)) ' ' num2str(tpsnr.cr(index(i)+d)) ']']);
            end;
            
        case 'B'

            % Get picture
            imB = RGBtoYCbCrimage( ...
                readsequence(seqpath, index(i) + d, index(i) + d), ...
                ss, 'bilinear');

            % BACKWARD PREDICTION FROM NEAREST I OR P PICTURE
            
            % Estimate motion vertical and horizontal vectors to reference
            [BKmvv BKmhv BKerrors] = estimatemotion(applypadding(imB.Y, 16), ...
                backwardimr.Y, [16 16], [15 15], 'fastfullsearch');

            % Apply motion estimated to the reference image with and
            BKimPestimated = applymotion(backwardimr.Y, BKmvv, BKmhv, [16 16]);
       
            % FORWARD PREDICTION FROM NEAREST I OR P PICTURE
            
            % Estimate motion vertical and horizontal vectors to reference
            [FRmvv FRmhv FRerrors] = estimatemotion(applypadding(imB.Y, 16), ...
                forwardimr.Y, [16 16], [15 15], 'fastfullsearch');

            % Apply motion estimated to the reference image with and
            FRimPestimated = applymotion(forwardimr.Y, FRmvv, FRmhv, [16 16]);
       
            % BIDERECTIONAL PREDICTION FROM NEARESTS I OR P PICTURES
            BIDimPestimated = uint8(0.5 * (double(BKimPestimated) + double(FRimPestimated)));
            BIDerrors = blkproc((double(BIDimPestimated) - double(applypadding(imB.Y, 16))) .^ 2, [16 16], 'sum(x(:))');
            
            % Calculate minimum prediction errors.
            predictionerrors = min(BKerrors, min(FRerrors, BIDerrors));

            % Mark what kind of prediction is used in each macroblock
            pmarker = zeros(size(predictionerrors));
            pmarker(BKerrors == predictionerrors) = 1;
            pmarker(FRerrors == predictionerrors) = 3;
            pmarker(BIDerrors == predictionerrors) = 5;

            % Kronecker product to resize pmarker.
            rpmarker = kron(pmarker, ones(16));
            
            % Get reference image and motion vectores.
            imBestimated = BKimPestimated;
            imBestimated(rpmarker == 3) = ...
                FRimPestimated(rpmarker == 3);
            imBestimated(rpmarker == 5) = ...
                BIDimPestimated(rpmarker == 5);
            
            % Differences between imB estimated and the original imB.
            imBdiff = double(applypadding(imB.Y, 16)) - ...
                double(imBestimated);
            
            % Encode diferences and motion vectores and prediction mask.
            diffcm = encodeimage(imBdiff, step2, 'motiondifference');
            maskBK = pmarker <= 2;
            mvcm = encodemv(BKmvv, BKmhv, 'motion2', ...
                FRmvv, FRmhv, maskBK | pmarker >= 5, ~maskBK);
            pmarkercm = encode(pmarker(:)', pmarkercode);
            
            Bcodedmessage = strcat(pmarkercm, mvcm, diffcm);
            
            % Write coded message in the payload segment.
            len = length(Bcodedmessage);
            while(len > cmsize - cmfilled)
                codedframesmessage(1, end + 1 : end + cmsize) = '0';
                cmsize = 2 * cmsize;
            end
            codedframesmessage(1, cmfilled + 1 : cmfilled + len) = ...
                Bcodedmessage(1, :);
            cmfilled = cmfilled + len;
            
            local_imBdiff = decoderesult(imBdiff, step2, 'difference');
            local_forwardimr.Y = uint8(double(imBestimated) + local_imBdiff);
            
            for j = 2 : 3
                
                backwardimrf = imfilter(...
                    getfield(backwardimr, fnames{j}), ...
                    h2, 'symmetric');
                BKimPestimated = applymotion( ...
                    getfield(backwardimr, fnames{j}), ...
                    BKmvv/2, BKmhv/2, [8 8]);
                
                FRimPestimated = applymotion( ...
                    getfield(forwardimr, fnames{j}), ...
                    FRmvv/2, FRmhv/2, [8 8]);
                
                BIDimPestimated = uint8(...
                    0.5 * (double(BKimPestimated) + double(FRimPestimated)));           
                
                % Kronecker product to resize pmarker.
                rpmarker = kron(pmarker, ones(8));
                
                % Get reference image and motion vectores.
                imBestimated = BKimPestimated;
                imBestimated(rpmarker == 3) = ...
                    FRimPestimated(rpmarker == 3);
                imBestimated(rpmarker == 5) = ...
                    BIDimPestimated(rpmarker == 5);
                
                imBdiff = ...
                    double(applypadding(getfield(imB, fnames{j}), 8)) - ...
                    double(imBestimated);
                
                diffcm = encodeimage(imBdiff, step2, 'motiondifference');
                
                len = length(diffcm);
                while(len > cmsize - cmfilled)
                    codedframesmessage(1, end + 1 : end + cmsize) = '0';
                    cmsize = 2 * cmsize;
                end
                codedframesmessage(1, cmfilled + 1 : cmfilled + len) = ...
                    diffcm(1, :);
                cmfilled = cmfilled + len;

                local_imBdiff = decoderesult(imBdiff, step2, 'difference');
                local_forwardimr = setfield(local_forwardimr, fnames{j}, ...
                    uint8(double(imBestimated) + local_imBdiff));
            
            end
     
            tnbits(index(i)+d) = cmfilled - cmfilled_old;
            ttype(index(i)+d) = 'B';
            tpsnr.y(index(i)+d) = mypsnr(local_forwardimr.Y,imB.Y);
            tpsnr.cb(index(i)+d) = mypsnr(local_forwardimr.Cb,imB.Cb);
            tpsnr.cr(index(i)+d) = mypsnr(local_forwardimr.Cr,imB.Cr);
            cmfilled_old = cmfilled;            
            if (verbose == 1),
                disp(['Coding ' ttype(index(i)+d) ' frame (' num2str(index(i)+d),'):' ...
                    ' bits=' num2str(tnbits(index(i)+d)) ...
                    ', psnr=[' num2str(tpsnr.y(index(i)+d)) ' ' num2str(tpsnr.cb(index(i)+d)) ' ' num2str(tpsnr.cr(index(i)+d)) ']']);
            end;

        otherwise
            error('Something went really wrong on picking picture type!!');
    end
end

% PACKBITS! AND CUT OF UNEEDED CODEFRAMESMESSAGE SPACE!!!
codedframes = packbits(codedframesmessage(1, 1 : cmfilled));
save(codedseqpath, 'codedframes');

if (verbose == 1),
    disp(strcat('TOTAL bits = ',num2str(cmfilled))); 
end;


function im = decoderesult(im, step, varargin)

if nargin > 2
    type = varargin{1};
else
    type = 'image';
end

switch type
    case 'image'
        im = blkproc(im, [8 8], @zigzagdctq8, step);
    case 'difference'
        im = blkproc(im, [8 8], @zigzagdctq8diff, step);
    otherwise
        error('Wrong type! Type may be: image or difference!');    
end

function im = applypadding(im, max)
padding = mod(max - mod(size(im), max), max);
im = im([1 : end end - 1 : -1 : end - padding(1, 1)], ...
    [1 : end end - 1 : -1 : end - padding(1, 2)]);

function im = zigzagdctq8(im, step)
im = uint8(idct2(iquantizedct8(quantizedct8(dct2(im), step), step)));

function imdiff = zigzagdctq8diff(imdiff, step)
imdiff = round(idct2(iquantizedct8diff(quantizedct8diff(dct2(imdiff), step), step)));