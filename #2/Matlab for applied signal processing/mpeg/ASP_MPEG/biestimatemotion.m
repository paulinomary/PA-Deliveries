function [imBestimated,imBdiff,pmarker] = biestimatemotion(im,prev_ref,past_ref,bs,sa,algorithm)
     
% BACKWARD PREDICTION
% Estimate motion vertical and horizontal vectors to reference image
[BKmvv BKmhv BKerrors] = estimatemotion(im, past_ref, bs, sa,algorithm);

% Apply motion estimated to the reference image
BKimPestimated = applymotion(past_ref, BKmvv, BKmhv, bs);

% FORWARD PREDICTION

% Estimate motion vertical and horizontal vectors to reference
[FRmvv FRmhv FRerrors] = estimatemotion(im, prev_ref, bs, sa,algorithm);

% Apply motion estimated to the reference image 
FRimPestimated = applymotion(prev_ref, FRmvv, FRmhv, bs);

% BIDERECTIONAL PREDICTION FROM NEARESTS I OR P PICTURES
BIDimPestimated = uint8(0.5 * (double(BKimPestimated) + double(FRimPestimated)));
switch(algorithm)
    case {'fullsearch','nstep'}
        BIDerrors = blkproc(abs(double(BIDimPestimated) - double(im)), bs, 'sum(x(:))');
    case 'fastfullsearch'
        BIDerrors = blkproc((double(BIDimPestimated) - double(im)) .^ 2, bs, 'sum(x(:))');
end

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
imBestimated(rpmarker == 3) = FRimPestimated(rpmarker == 3);
imBestimated(rpmarker == 5) = BIDimPestimated(rpmarker == 5);

% Differences between imB estimated and the original imB.
imBdiff = double(im) - double(imBestimated);

