function h = encodeBchoices(seqpath, first, last, varargin)

% Set the period of P and I pictures.
pP = 5;
pI = 4;
if nargin > 3
    pP = varargin{1};
    if nargin > 4
        pI = varargin{2};
    end
end

% Set number of frames.
nframes = last - first + 1;

% Set deviation from 1 to first.
d = first - 1;

% Get sequence size.
seqsize = sequencesize(seqpath);

% Set stream codification order.
[pictype index] = shufflestream(pP, pI, nframes);

% Factor used to operate the diferent Y' and Cb', Cr' picture componentes.
cf = [16 8 8];

% Sub-sampling factor.
ss = [2 2];

% Steps used in DCT quantization.
step1 = jpegsteps;
step2 = jpegsteps;
step2(1, 1) = step2(1, 1) * 2;

% Low-Pass filters
h = [1 2 1; 2 4 2; 1 2 1] / 16;
h2 = [0 1 0; 1 4 1; 0 1 0] / 8;

% Y' Cb' Cr' structure.
YCbCrstruct = struct('Y', [], 'Cb', [], 'Cr', [], 'ss', ss);
fnames = fieldnames(YCbCrstruct);

% Backward and foward reference image.
backwardimr = YCbCrstruct;
forwardimr = YCbCrstruct;

h(pP, pI, 1 : 6) = 0;

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
                % Update forward image reference with decoder results.
                forwardimr = setfield(forwardimr, fnames{j}, ...
                    decoderesult(applypadding( ...
                    getfield(imI, fnames{j}), cf(1 ,j)), step1));
            end
            
        case 'P'

            % Update backward image reference.
            backwardimr = forwardimr;
                        
            % Get picture
            imP = RGBtoYCbCrimage( ...
                readsequence(seqpath, index(i) + d, index(i) + d), ...
                ss, 'bilinear');

            % Apply low-pass filter to reference image.
            backwardimrfY = imfilter(backwardimr.Y, h, 'symmetric');
            
            % Estimate motion vertical and horizontal vectors to reference
            % image with and without filter.
            [mvvf mhvf errorsf] = estimatemotion(applypadding(imP.Y, 16), ...
                backwardimrfY, [16 16], [15 15], 'fastfullsearch');
            [mvv mhv errors] = estimatemotion(applypadding(imP.Y, 16), ...
                backwardimr.Y, [16 16], [15 15], 'fastfullsearch');

            % Apply motion estimated to the reference image with and
            % without filter.
            imPestimatedf = applymotion(backwardimrfY, mvvf, mhvf, [16 16]);
            imPestimated = applymotion(backwardimr.Y, mvv, mhv, [16 16]);
       
            % Predicion mask decides wether to use filtered predicion or
            % not.
            pmask = errorsf < errors;
            
            % Kronecker product to resize prediction mask.
            rpmask = logical(kron(pmask, ones(16)));
        
            % Get diferences used and motion vectores.
            imPestimated(rpmask) = imPestimatedf(rpmask);
            mvv(pmask) = mvvf(pmask);
            mhv(pmask) = mhvf(pmask);

            % Differences between imP estimated and the original imP.
            imPdiff = double(applypadding(imP.Y, 16)) - ...
                double(imPestimated);
            
            % Apply decoder quantization errors to differences.
            imPdiff = decoderesult(imPdiff, step2, 'difference');
            
            % Sum the differences to the image estimated with motion
            % vectors.
            forwardimr.Y = double(imPestimated) + imPdiff;
            
            % Compute and encode Cb' and Cr' components differences.
            for j = 2 : 3

                backwardimrf = imfilter(...
                    getfield(backwardimr, fnames{j}), ...
                    h2, 'symmetric');
                
                imPestimatedf = applymotion(... 
                    backwardimrf, mvvf/2, mhvf/2, [8 8]);
                imPestimated = applymotion( ...
                    getfield(backwardimr, fnames{j}), mvv/2, mhv/2, [8 8]);

                rpmask = logical(kron(pmask, ones(8)));
                imPestimated(rpmask) = imPestimatedf(rpmask);
                        
                imPdiff = ...
                    double(applypadding(getfield(imP, fnames{j}), 8)) - ...
                    double(imPestimated);
            
                imPdiff = decoderesult(imPdiff, step2, 'difference');
                forwardimr = setfield(forwardimr, fnames{j}, ...
                    double(imPestimated) + imPdiff);
            end
            
        case 'B'
        
            % Get picture
            imB = RGBtoYCbCrimage( ...
                readsequence(seqpath, index(i) + d, index(i) + d), ...
                ss, 'bilinear');

            % BACKWARD PREDICTION FROM NEAREST I OR P PICTURE
            % Apply low-pass filter to reference image.
            backwardimrfY = imfilter(backwardimr.Y, h, 'symmetric');
            
            % Estimate motion vertical and horizontal vectors to reference
            % image with and without filter.
            [BKmvvf BKmhvf BKerrorsf] = estimatemotion(applypadding(imB.Y, 16), ...
                backwardimrfY, [16 16], [15 15], 'fastfullsearch');
            [BKmvv BKmhv BKerrors] = estimatemotion(applypadding(imB.Y, 16), ...
                backwardimr.Y, [16 16], [15 15], 'fastfullsearch');

            % Apply motion estimated to the reference image with and
            % without filter.
            BKimPestimatedf = applymotion(backwardimrfY, BKmvvf, BKmhvf, [16 16]);
            BKimPestimated = applymotion(backwardimr.Y, BKmvv, BKmhv, [16 16]);
       
            % FORWARD PREDICTION FROM NEAREST I OR P PICTURE
            % Apply low-pass filter to reference image.
            forwardimrfY = imfilter(forwardimr.Y, h, 'symmetric');
            
            % Estimate motion vertical and horizontal vectors to reference
            % image with and without filter.
            [FRmvvf FRmhvf FRerrorsf] = estimatemotion(applypadding(imB.Y, 16), ...
                forwardimrfY, [16 16], [15 15], 'fastfullsearch');
            [FRmvv FRmhv FRerrors] = estimatemotion(applypadding(imB.Y, 16), ...
                forwardimr.Y, [16 16], [15 15], 'fastfullsearch');

            % Apply motion estimated to the reference image with and
            % without filter.
            FRimPestimatedf = applymotion(forwardimrfY, FRmvvf, FRmhvf, [16 16]);
            FRimPestimated = applymotion(forwardimr.Y, FRmvv, FRmhv, [16 16]);
       
            % BIDERECTIONAL PREDICTION FROM NEARESTS I OR P PICTURES
            BIDimPestimatedf = 0.5 * (BKimPestimatedf + FRimPestimatedf);
            BIDimPestimated = 0.5 * (BKimPestimated + FRimPestimated);           
            BIDerrorsf = blkproc((BIDimPestimatedf - double(applypadding(imB.Y, 16))) .^ 2, [16 16], 'sum(x(:))');
            BIDerrors = blkproc((BIDimPestimated - double(applypadding(imB.Y, 16))) .^ 2, [16 16], 'sum(x(:))');
            
            predictionerrors = min(BKerrorsf, ...
                min(BKerrors, min(FRerrorsf, ...
                min(FRerrors, min(BIDerrorsf, BIDerrors)))));

            p = ones(size(predictionerrors)) * -1;
            p(BKerrors == predictionerrors) = 1;
            p(BKerrorsf == predictionerrors) = 2;
            p(FRerrors == predictionerrors) = 3;
            p(FRerrorsf == predictionerrors) = 4;
            p(BIDerrors == predictionerrors) = 5;
            p(BIDerrorsf == predictionerrors) = 6;

            h(pP, pI, 1) = h(pP, pI, 1) + length(p(p == 1));
            h(pP, pI, 2) = h(pP, pI, 2) + length(p(p == 2));
            h(pP, pI, 3) = h(pP, pI, 3) + length(p(p == 3));
            h(pP, pI, 4) = h(pP, pI, 4) + length(p(p == 4));
            h(pP, pI, 5) = h(pP, pI, 5) + length(p(p == 5));
            h(pP, pI, 6) = h(pP, pI, 6) + length(p(p == 6));
             

        otherwise
            error('Something went really wrong on picking picture type!!');
    end
end

function im = decoderesult(im, step, varargin)

if nargin > 2
    type = varargin{1};
else
    type = 'image';
end

switch type
    case 'image'
        im = double(blkproc(im, [8 8], @zigzagdctq8, step));
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
im = idct2(iquantizedct8(quantizedct8(dct2(im), step), step));

function imdiff = zigzagdctq8diff(imdiff, step)
imdiff = idct2(iquantizedct8(quantizedct8(dct2(imdiff), step), step));

