function mpegdecode(codedseqpath, lQk, peI, peP, peB)

Qk = 1.0;
if nargin > 1,
    Qk = lQk;
end

% Load coded segment and unpackbits.
load(codedseqpath);
codedframesmessage = unpackbits(codedframes);

% Extract header from the codedframesmessage, and its information:
% Get sequence size
seqsize = [bin2dec(codedframesmessage(1, 1 : 9)), ...
        bin2dec(codedframesmessage(1, 10 : 18))];
% Get number of frames to decode
nframes = bin2dec(codedframesmessage(1, 19 : 27));
% Get period of P and I pictures
pP = bin2dec(codedframesmessage(1, 28 : 34));
pI = bin2dec(codedframesmessage(1, 35 : 40));

dlen = 500000;

if dlen > length(codedframesmessage)
    codedframesmessage = codedframesmessage(1, 40 + 1 : end);
    lefttodecode = '';
    dlen = 0;
else
    lefttodecode = codedframesmessage(1, 40 + dlen + 1 : end);
    codedframesmessage = codedframesmessage(1, 40 + 1 : 40 + dlen);
end

% Set stream codification order.
[pictype index] = shufflestream(pP, pI, nframes);

% Sub-sampling factor
ss = [2 2];

% Load decoder used to decode predictive mask.
load('./decoders/pcode/pdecoder.mat');

% Steps used in DCT quantization.
step1 = jpegsteps.*Qk;
%step2 = jpegsteps.*Qk; step2(1, 1) = step2(1, 1) * 2;
step2 = 16*ones(size(step1)).*Qk*2;

% Padded image size for Y' component.
pimsize = seqsize + mod(16 - mod(seqsize, 16), 16);

% Predictive mask length.
pmasklen = prod(pimsize / 16);

% Low-Pass filters
h = [1 2 1; 2 4 2; 1 2 1] / 16;
h2 = [0 1 0; 1 4 1; 0 1 0] / 8;

% Y' Cb' Cr' structure.
YCbCrstruct = struct('Y', [], 'Cb', [], 'Cr', [], 'ss', ss);
fnames = fieldnames(YCbCrstruct);

% Backward and foward reference image.
backwardimr = YCbCrstruct;
forwardimr = YCbCrstruct;
imtowrite = YCbCrstruct;

for i = 1 : nframes

    if dlen > length(codedframesmessage)
        if length(lefttodecode) < dlen
            dlen = length(lefttodecode);
        end
        codedframesmessage(1, end + 1 : end + dlen) = lefttodecode(1, 1 : dlen);
        lefttodecode = lefttodecode(1, dlen + 1 : end);
    end
        
    switch pictype(index(i))
        case 'I'
            % Update backward image reference.
            backwardimr = forwardimr;
            
            % Decode image.
            [imI, codedframesmessage] = ...
                decodeimage(codedframesmessage, pimsize, step1);
            
            % Update forward image reference and image that is going to be
            % writen.
            forwardimr.Y = imI;
   %        imtowrite.Y = removepadding(uint8(imI), seqsize);
            
            for j = 2 : 3
                [imI, codedframesmessage] = ...
                    decodeimage(codedframesmessage, pimsize / 2, step1);
                forwardimr = setfield(forwardimr, fnames{j}, imI);
   %            imtowrite = setfield(imtowrite, fnames{j}, ...
   %               removepadding(uint8(imI), ceil(seqsize / 2)));
            end
    
            if peI > 0
               forwardimr = addnoise(forwardimr, peI); 
            end
            
            imtowrite.Y = removepadding(uint8(forwardimr.Y), seqsize);
            imtowrite.Cb = removepadding(uint8(forwardimr.Cb), ceil(seqsize/2));
            imtowrite.Cr = removepadding(uint8(forwardimr.Cr), ceil(seqsize/2));
            
            
        case 'P'
            % Update backward image reference.
            backwardimr = forwardimr;
                                    
            % Decode motion vectors
            [mvv, mhv, codedframesmessage] = ...
                decodemv(codedframesmessage(1, 1 : end), ...
                pimsize, [16 16]);
                        
            % Apply motion vectors to the reference image with and
            imPestimated = applymotion(backwardimr.Y, mvv, mhv, [16 16]);
            
            % Decode differences between the esimated and original images.
            [imPdiff codedframesmessage] = decodeimage( ...
                codedframesmessage, pimsize, step2, 'motiondifference');
            
            % Sum the differences to predictive image.
            imP = uint8(double(imPestimated) + imPdiff);
            
            
            % Update forwardimr and save imP to imtowrite.
            forwardimr.Y = imP;
            %imtowrite.Y = removepadding(uint8(imP), seqsize);
            
            for j = 2 : 3
                imPestimated = applymotion( ...
                    getfield(backwardimr, fnames{j}), ...
                    mvv/2, mhv/2, [8 8]);
                
                [imPdiff codedframesmessage] = decodeimage( ...
                    codedframesmessage, pimsize / 2, step2, ...
                    'motiondifference');
                imP = uint8(double(imPestimated) + imPdiff);
                
                forwardimr = ...
                    setfield(forwardimr, fnames{j}, imP);
             %   imtowrite = setfield(imtowrite, fnames{j}, ...
              %      removepadding(uint8(imP), ceil(seqsize / 2)));
            end
            
            if peP > 0
               forwardimr = addnoise(forwardimr, peP); 
            end
            
            imtowrite.Y = removepadding(uint8(forwardimr.Y), seqsize);
            imtowrite.Cb = removepadding(uint8(forwardimr.Cb), ceil(seqsize/2));
            imtowrite.Cr = removepadding(uint8(forwardimr.Cr), ceil(seqsize/2));
            
            
            
        case 'B'
            % Decode prediction marker.
            [pmarker codedframesmessage] = ...
                decodemessage(codedframesmessage, pdecoder, pmasklen);
            pmarker = reshape(pmarker, pimsize / 16);
            
            % Decode motion vectors.
            maskBK = pmarker <= 2;
            [BKmvv BKmhv FRmvv FRmhv codedframesmessage] = decodemv( ...
                codedframesmessage, pimsize, [16 16], 'motion2', ...
                maskBK | pmarker >= 5, ~maskBK);
            
            % BACKWARD PREDICTION FROM NEAREST I OR P PICTURE
            
            % Apply motion estimated to the reference image with and
            % without filter.
            BKimPestimated = double( ...
                applymotion(backwardimr.Y, BKmvv, BKmhv, [16 16]));
       
            % FORWARD PREDICTION FROM NEAREST I OR P PICTURE
            
            % Apply motion estimated to the reference image with and
            FRimPestimated = double( ...
                applymotion(forwardimr.Y, FRmvv, FRmhv, [16 16]));
       
            % BIDERECTIONAL PREDICTION FROM NEARESTS I OR P PICTURES
            BIDimPestimated = 0.5 * (BKimPestimated + FRimPestimated);        
                
            % Kronecker product to resize prediction marker.
            rpmarker = kron(pmarker, ones(16));
            
            % Get reference image and motion vectores.
            imBestimated = BKimPestimated;
            imBestimated(rpmarker == 3) = ...
                FRimPestimated(rpmarker == 3);
            imBestimated(rpmarker == 5) = ...
                BIDimPestimated(rpmarker == 5);
            
            % Decode differences between the esimated and original images.
            [imBdiff codedframesmessage] = decodeimage( ...
                codedframesmessage, pimsize, step2, 'motiondifference');
            
            % Sum the differences to predictive image.
            imB = uint8(double(imBestimated) + imBdiff);
            
            % Update forwardimr and save imP to imtowrite.
            %imtowrite.Y = removepadding(uint8(imB), seqsize);
            imtowrite.Y = imB;
            
            for j = 2 : 3
                BKimPestimated = double( ...
                    applymotion(getfield(backwardimr, fnames{j}), ...
                    BKmvv/2, BKmhv/2, [8 8]));
                
                FRimPestimated = double( ...
                    applymotion(getfield(forwardimr, fnames{j}), ...
                    FRmvv/2, FRmhv/2, [8 8]));
                
                BIDimPestimated = ...
                    0.5 * (BKimPestimated + FRimPestimated);        
                
                rpmarker = kron(pmarker, ones(8));
                imBestimated = BKimPestimated;
                imBestimated(rpmarker == 3) = ...
                    FRimPestimated(rpmarker == 3);
                imBestimated(rpmarker == 5) = ...
                    BIDimPestimated(rpmarker == 5);
                
                [imBdiff codedframesmessage] = decodeimage( ...
                    codedframesmessage, pimsize / 2, step2, ...
                    'motiondifference');
                
                imB = uint8(double(imBestimated) + imBdiff);
                
                %imtowrite = setfield(imtowrite, fnames{j}, ...
                    %removepadding(uint8(imB), ceil(seqsize / 2)));

                imtowrite = setfield(imtowrite, fnames{j}, imB);
            end
            
            if peB > 0
               imtowrite = addnoise(imtowrite, peB); 
            end
            
            imtowrite.Y = removepadding(uint8(imtowrite.Y), seqsize);
            imtowrite.Cb = removepadding(uint8(imtowrite.Cb), ceil(seqsize/2));
            imtowrite.Cr = removepadding(uint8(imtowrite.Cr), ceil(seqsize/2));
            
        otherwise
            error('Something went really wrong on picking picture type!!');
    end
    
    % Write image to disk.
    if index(i) - 1 < 10
        imwrite(YCbCrtoRGBimage(imtowrite, 'bilinear'), ...
            strcat(codedseqpath, '00', num2str(index(i) - 1), '.bmp'), 'BMP');
    elseif index(i) - 1 < 100
        imwrite(YCbCrtoRGBimage(imtowrite, 'bilinear'), ...
            strcat(codedseqpath, '0', num2str(index(i) - 1), '.bmp'), 'BMP');
    else
        imwrite(YCbCrtoRGBimage(imtowrite, 'bilinear'), ...
            strcat(codedseqpath, num2str(index(i) - 1), '.bmp'), 'BMP');
    end
end

function im = removepadding(im, imsize)
im = im(1 : imsize(1), 1 : imsize(2));


%%%%%%%%%%%%%%%%%%% NEW %%%%%%%%%%%%%%
function imYCbCr = addnoise(imYCbCr, pe)
    blkerror = rand(size(imYCbCr.Y)/8) < pe;
    blkerror8 = kron(blkerror, ones(8,8));
    blkerror4 = kron(blkerror, ones(4,4));

    errorY = 40;
    errorCb = 10;
    errorCr = 20;
    
    imYCbCr.Y(logical(blkerror8)) = imYCbCr.Y(logical(blkerror8)) + errorY;
    imYCbCr.Cr(logical(blkerror4)) = imYCbCr.Cb(logical(blkerror4)) + errorCb;
    imYCbCr.Cr(logical(blkerror4)) = imYCbCr.Cr(logical(blkerror4)) + errorCr;    
    