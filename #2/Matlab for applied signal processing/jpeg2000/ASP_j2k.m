

%% Chapter XX - How does digital cinema compresses images ?
%
% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%
% This file uses the Wavelet Toolbox of MATLAB.

clear all
close all
nfig=1;
save_fig=0; % '1' if figures have to be saved in fig and pdf format in the 'figures' folder; '0' otherwise

%% 1. Experiments with the Wavelet Transform
% In this section four numerical experiments are presented to make the
% reader more familiar with the practical aspects of the 1-D and 2-D DWT
% implementations in Matlab$^\copyright$. The third experiment goes
% slightly further this simple exploration by realizing a very simple
% image compression procedure in the wavelet domain.

%% 1.1 Experiment 1 : Computing Haar DWT
% Let us realize an example of DWT computation thanks to the use of the
% Matlab Wavelet Toolbox$^{\rm TM}$ \cite{matlabwavtbx}. Haar wavelet
% transform is there very simple to use. The following sequence of Matlab
% commands performs a DWT of our favorite 1-D signal (Fig.
% \ref{fig:ex-sig-1-d}) from resolution $J=10$ to resolution $J_0=5$.
% Detail and approximation coefficients are represented with the {\tt
% waveshow()} command. Its output is given in the following figure.
% Localization of important coefficients close to transient signal parts is
% even more clear now than in the previous example.

%Loading the signal, N=1024=2^10
load 1d-sig;
% Performing the DWT,
J = 10; J0 = 5;
[W,L] = wavedec(sig, J-J0, 'haar');
% Showing it (waveshow.m is part of the CDROM)
h=figure; waveshow(sig,W,L)

if save_fig==1
    saveas(h,sprintf('figures/fig%d.fig',nfig));
    saveas(h,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end

% Showing the ladder of approximations (appshow.m is part of the CDROM)
h=figure; appshow(sig,W,L, 'haar')

if save_fig==1
    saveas(h,sprintf('figures/fig%d.fig',nfig));
    saveas(h,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end

%% 1.2 Experiment 2 : Performing 2D-DWT
% The following routine computes the 2-D DWT of the {\tt Cameraman} picture provided in
% Matlab. The code is given below. It displays the wavelet coefficients in the same
% shape than the one of the the matrix ${\bf W}$ above. Horizontal,
% vertical and diagonal coefficients, which are normalized in {\tt
% waveshow2()} by their maximal absolute value at each resolution,
% display the same kind of behavior than for 1-D signal decomposition :
% there are important in amplitude only in the transient parts,
% i.e. close to the edges or to the textures of the image (e.g. the
% grass in the bottom).

% Loading the image
% img is a 256x256 size array 
%im = double(imread('cameraman.tif'));
im = double(imread('images/Barbara256.pgm'));
h=figure; imshow(im,[7 253]);

if save_fig==1
    saveas(h,sprintf('figures/fig%d.fig',nfig));
    saveas(h,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end

% 2D DWT Computations with Daubechies 9/7 (== 'bior4.4')
% W contains the DWT coefficients in a column shape.
J = log2(256); J0 = 2; wname = 'bior4.4';
dwtmode('per');
[W,S] = wavedec2(im, J-J0, wname);
% The DWT array
h=figure; waveshow2(W,S,wname);

if save_fig==1
    saveas(h,sprintf('figures/fig%d.fig',nfig));
    saveas(h,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end


%% 1.3 Experiment 3 : On the compression road.
% As a new numerical example, we  want to insist now on the 
% ``sparseness'' of the Wavelet coefficients for the representation 
% of \emph{natural} images, i.e. made of structured content as objects 
% separated by edges. 
% 
% The \emph{compression principle} of the WT holds in this : \textbf{a very
% few  number of coefficients concentrate the essential of the image
% information}. This can already be perceived the previous Figure 
% \ref{fig:2d-dwt-example} [[FIX THE FIGURE LABEL!]]. Conversely to initial
% pixel values, detail coefficients with high amplitude are not very
% numerous and well localized on image edges.
% 
% The following code is the continuation of the code from Expermient 2. 
% The idea is to realize a compression of
% 90\% of the {\tt cameraman} image by keeping only 10\% of its
% strongest wavelet coefficients. This is achieved by sorting the
% wavelet values by decreasing order of magnitude (L3) and recording the
% amplitude of the $K$th strongest element (with $K$ the closest integer
% to $N^2/10$). Then, all the wavelet coefficients with magnitude are
% set to~0 (L10--11), i.e. their information is lost. The rebuilding of
% the image is then realized with the new values in L8. The result is
% presented in the following figure.

% Sorting wavelet coefficient amplitude
sW = sort(abs(W(:)), 'descend');

% Number of elements to keep
% Compression of 90% !
K = round(256*256/10);
T = sW(K);

% Thresholding values of W lesser than T
% i.e. we keep the K strongest
nW = W;
nW(abs(nW) < T) = 0;

% Rebuilding the compressed image
Timg = waverec2(nW, S, wname);
h=figure; imshow(Timg,[7 253]);

if save_fig==1
    saveas(h,sprintf('figures/fig%d.fig',nfig));
    saveas(h,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end

%% 1.4 Experiment 4 : Quantifying Compression Quality
% Let us now quantify a bit further the quality reached by the compression
% scheme of Experiment 3 in function of both the number of DWT coefficients
% kept during the thresholding and the type of wavelet filters used. 
% 
% In the following figure 
% [[Figure: quality_compr.pdf, 
%   Caption: Quality curve of compressed
%            images (for the {\tt Cameraman} image) for different
%            percentage of DWT coefficients and for different filters. 
%   Label: fig:qual-curve]]
% the quality curve obtained for different percentage of wavelet
% coefficients is drawn (from 5\% to 30\%, i.e. from 95\% to 70\% of
% compression) for the {\tt Cameraman} image of Figure
% \ref{fig:camer-orig} [[FIX THE FIGURE LABEL!]]. We can clearly see that the Daubechies 9/7 filter
% provides the best quality. However, quality is not the only criterion
% which makes a filter better than another. Daubechies (or Legall) 5/3
% filters, which can be expressed by rational numbers, are used for
% \emph{Lossless} compression in JPEG 2000. The irrationality of the
% Daubechies 9/7 filters prevents this possibility. <comment: JPEG-2000
% experts have to add something here ;-)> <comment: code below is perhaps
% not really important>


nbpix = 256*256;

% Fixing the percentage of pixels to keep in the compression
% between 5% and 30%
K = round(((5:5:30)/100) * nbpix);
nbK = length(K);

% Image decomposed on J levels
%im = double(imread('cameraman.tif'));     
im = double(imread('images/Barbara256.pgm'));
J = log2(256); J0 = 0; 

% Quality metrics between two images : MSE and PSNR.
% Assuming an 8 bits original image
MSE = @(X,Y) norm(X(:) - Y(:), 2)^2 / nbpix; 
PSNR = @(X,Y) 10*log10( (256-1)^2 / MSE(X,Y) );

wname1 = 'bior4.4'; %% Daubechies 9/7
wname2 = 'bior3.3'; %% Daubechies/Legall 5/3

[W1,S1] = wavedec2(im, J-J0, wname1);
[W2,S2] = wavedec2(im, J-J0, wname2);

sW1 = sort(abs(W1(:)), 'descend');
sW2 = sort(abs(W2(:)), 'descend');

for k = 1:nbK,
   % Setting all the DWT coefficients smaller than the Kth magnitude to
   % zero.
   % For DB97
   T1 = sW1(K(k));
   nW1 = W1;
   nW1(abs(nW1) < T1) = 0;
   Timg1 = waverec2(nW1, S1, wname1);

   % For DB53
   T2 = sW2(K(k));
   nW2 = W2;
   nW2(abs(nW2) < T2) = 0;
   Timg2 = waverec2(nW2, S2, wname2);

   % Recording quality
   curve_97(k) = PSNR(im, Timg1); 
   curve_53(k) = PSNR(im, Timg2);     
end

% Plotting the result
h=figure;
plot(100*K/nbpix, curve_97, 'o-', 100*K/nbpix, curve_53, 'o--','LineWidth',1);
axis tight
title('PSNR versus number of DWT coeffcients');
xlabel('Percentage of coefficients (w.r.t. the number of initial pixels).');
ylabel('Quality (in dB)');
legend('Daub. 9/7', 'Daub. 5/3', 'Location', 'SouthEast');
grid on;

if save_fig==1
    saveas(h,sprintf('figures/fig%d.fig',nfig));
    saveas(h,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end

%% 2. A simplified JPEG 2000 scheme
% In this second part of the code, an image is chosen and processed through
% a simplified JPEG 2000 coding scheme. Some image examples are given but
% any other grayscale image could be specified (color images could also of
% course be supported but require minor changes in the code).

%% Parameters
% These parameters will be used through the whole simplified JPEG 2000
% coding scheme. Practically, a structure IMG is defined that will contain
% the required settings to process the image. A small description is given
% beside each field of the structure. The Context Distribution Table will
% be further explained below.

%img.path = 'images/lena512.pgm';       % path to reach the image to be processed.
img.path = 'images/Barbara512.pgm';
%img.path = 'cameraman.tif';
img.useownCDT = 1;                      % '1' if the CDT computed on the processed image should be used ; '0' otherwise
img.CDT2use = 'images/Barbara512_CDT.mat';   % path to the Context Distribution Table to be used if img.useownCDT=0;
img.wfilt = 0;                          % 0 = 5-3 wavelet transf., 1 = 9-7 wavelet transf.
img.nwdec = 5;                          % number of decompositions = number of res. levels - 1
img.cbh = 32;                           % code-block size
img.cbw = 32;

%% Pre-processing
% Before applying the DWT on the image, some pre-processing is achieved.
% This includes (1) an inter-component decorrelation (RGB => YCbCr) in case of
% color images, (2) the optional division of the image in tiles, and (3)
% the shifting of coefficients from an unsigned to a signed representation.
% In the following, we load a grayscale image (no inter-component
% decorrelation is performed), and no tile tiles are used, so that the
% pre-processing only consists in a DC-level shifting.

X = imread(img.path);
X = double(X);
img.bdepth = ceil(log2(max(X(:)+1)));
[img.h img.w] = size(X);
X = X - pow2(img.bdepth-1); % DC-level shifting

%% 2.1 Discrete Wavelet Transform
% The JPEG 2000 standard defines 2 different filters, namely the 5-3 and
% the 9-7 transform. The former is used in JPEG 2000 when lossless coding
% is required (to be accurate, the "integer-to-integer" version of this
% transform is used, so that it is truly reversible). The latter, with a
% slightly higher decorrelating power, is used for lossy coding.

%%
% We first define the two transforms, specifying the filter coefficients of
% the 5-3 and 9-7 tranformation, respectively. In the JPEG 2000 standard~\cite{normeJ2K}, a
% lifting scheme implementation is described so that the following
% coefficients are not explicitly mentioned. However, they can easily be
% found in~\cite{skodras2001jus} for example. The norm of the synthesis
% filters are also given. It will be used when approximating the distortion
% reduction brought in the pixel domain by the refinement of each wavelet
% coefficient. The two matrices wnorm_53 and wnorm_97 are structured as
% follows : rows correspond to each subband (1=LL, 2=HL, 3=LH, 4=HH) while
% columns correspond to successive resolution levels (1=highest resolution
% level, etc).

% 5-3 transformation
lo_53_D = [0 -1/8 2/8 6/8 2/8 -1/8];
hi_53_D = [0 -1/2 1 -1/2 0 0];
lo_53_R = [0 1/2 1 1/2 0 0];
hi_53_R = [0 -1/8 -2/8 6/8 -2/8 -1/8];

% norm of the 5-3 synthesis filter

wnorm_53 = [...
	1.000 1.500 2.750 5.375 10.68 21.34 42.67 85.33 170.7 341.3;...
	1.038 1.592 2.919 5.703 11.33 22.64 45.25 90.48 180.9 0;...
	1.038 1.592 2.919 5.703 11.33 22.64 45.25 90.48 180.9 0;...
	.7186 .9218 1.586 3.043 6.019 12.01 24.00 47.97 95.93 0 ...
    ];

% 9-7 transformation
lo_97_D = [0 0.02674875741080976 -0.01686411844287495 -0.07822326652898785 0.2668641184428723 0.6029490182363579 ...
    0.2668641184428723 -0.07822326652898785 -0.01686411844287495 0.02674875741080976];
hi_97_D = [0 0.09127176311424948 -0.05754352622849957 -0.5912717631142470 1.115087052456994 -0.5912717631142470 ...
    -0.05754352622849957 0.09127176311424948 0 0];
lo_97_R = [0 -0.09127176311424948 -0.05754352622849957 0.5912717631142470 1.115087052456994 0.5912717631142470 ...
    -0.05754352622849957 -0.09127176311424948 0 0];
hi_97_R = [0 0.02674875741080976 0.01686411844287495 -0.07822326652898785 -0.2668641184428723 0.6029490182363579 ...
    -0.2668641184428723 -0.07822326652898785 0.01686411844287495 0.02674875741080976];

% norm of the 9-7 synthesis filter

wnorm_97 = [...
	1.000 1.965 4.177 8.403 16.90 33.84 67.69 135.3 270.6 540.9;...
	2.022 3.989 8.355 17.04 34.27 68.63 137.3 274.6 549.0 0;...
	2.022 3.989 8.355 17.04 34.27 68.63 137.3 274.6 549.0 0;...
	2.080 3.865 8.307 17.18 34.71 69.59 139.3 278.6 557.2 0 ...
    ];

%%
% According to the parameter WFILT, we apply the 5-3 or 9-7 transform. The
% function wavedec2 from the Wavelet Toolbox is used to do this.

if img.wfilt==0
    [C,S] = wavedec2(X,img.nwdec,lo_53_D,hi_53_D);                  % dyadic decomposition
elseif img.wfilt==1
    [C,S] = wavedec2(X,img.nwdec,lo_97_D,hi_97_D);                  % dyadic decomposition
else
    error('wavelet filter not recognized');
end

img.S = S; % we record this matrix to re-use it when reconstructing the image

%% 2.2 Context-based modeling of coefficients bit-planes

% After the DWT, the JPEG 2000 algorithm performs a context-based entropy
% coding of each sub-band. This section illustrates the advantage of
% context-based modeling by comparing two measures of the entropy
% associated with the binary representation of the image DWT coefficients.
% In the first case, the binary symbols are assumed to be generated by an
% iid sequence of random variables, and the probability distribution of
% each binary random variable is estimated based on the frequency of
% occurence of one and zero symbols. In the second case, the probabilities
% of binary random variables are estimated conditionally to their context.
% The reduction of entropy between the first and second case corresponds to
% the benefit obtained from the chosen context model.

%%
% In the meantime, this section computes the incremental bit-budget and
% distortion reduction resulting from the addition of bitplanes to refine
% DWT coefficients. Thereby, it provides the inputs required by the
% rate-distortion optimal bit allocation mechanisms envisioned in the next
% section.

%%
% It should be noted that as we wanted to focus on the original aspects of
% JPEG 2000 and avoid a too long section, we did not actually implement the
% entropy coder, which compresses the input bit-stream according to the
% chosen probability distribution. No compressed image is therefore
% produced. Only the performance of such coder is evaluated (through
% estimation of the source entropy), such that the output rate for a given
% compression ratio can be estimated.


%% 2.2.1 Quantization and ranging
% Before being entropy-coded, wavelet coefficients are quantized and mapped
% on a certain amount of bits (ranging).

%%
% To do so, we first separate the sign and magnitude of the wavelet
% coefficients. They are indeed encoded separately in JPEG 2000.

Csign = sign(C);
Cmagn = abs(C);

%%
% In case of a lossless compression (implying an integer-to-integer 5-3
% transform), the quantization step is of course set to 1. In case of a
% lossy compression using the 9-7 transform, the quantization stepsize
% chosen in this example follows the rule used in the OpenJPEG
% library~\cite{openjpeg} : it makes it depend on the kind of subband
% involved and on the norm of the synthesis filter (which depends itself on
% the resolution level and on the subband).

%%
% When the coefficients have been quantized, they are mapped onto a certain
% amount of bits (fixed-point representation). In this simple experiment,
% we simply keep the integer part of the coefficients, that we represent
% using 16 bits, which is sufficient to avoid any overflow with images with
% reasonable bitdepths (between 8 and 12). In a true JPEG 2000 coding
% scheme, this number is of course computed based on the actual image
% bitdepth, so as to avoid any overflow due to the possible increase of the
% dynamic range following the DWT.

%%
% All quantized coefficients are stored in the IMG
% structure, resolution per resolution, and subband per subband.

% LL-subband
As = reshape(Csign(1:(S(1,1)*S(1,2))),S(1,1),S(1,2));
Am = reshape(Cmagn(1:(S(1,1)*S(1,2))),S(1,1),S(1,2));
if img.wfilt==0
    img.res(1).sb(1).qstep = 1;
elseif img.wfilt==1
    img.res(1).sb(1).qstep = 1/wnorm_97(1,img.nwdec+1);
end
img.res(1).sb(1).sign = As;
img.res(1).sb(1).coeff = uint16(floor(Am/img.res(1).sb(1).qstep));

% High frequencies subbands
for i=2:img.nwdec+1

    %coeff and sign extraction
    [Hm,Vm,Dm] = detcoef2('all',Cmagn,S,img.nwdec+2-i);
    [Hs,Vs,Ds] = detcoef2('all',Csign,S,img.nwdec+2-i);

    % HL subband
    if img.wfilt==0
        img.res(i).sb(1).qstep = 1;
    elseif img.wfilt==1
        img.res(i).sb(1).qstep = 2/wnorm_97(2,img.nwdec+2-i);
    end
    img.res(i).sb(1).sign = Hs;
    img.res(i).sb(1).coeff = uint16(floor(Hm/img.res(i).sb(1).qstep));

    % LH subband
    if img.wfilt==0
        img.res(i).sb(2).qstep = 1;
    elseif img.wfilt==1
        img.res(i).sb(2).qstep = 2/wnorm_97(3,img.nwdec+2-i);
    end
    img.res(i).sb(2).sign = Vs;
    img.res(i).sb(2).coeff = uint16(floor(Vm/img.res(i).sb(2).qstep));

    % HH subband
    if img.wfilt==0
        img.res(i).sb(3).qstep = 1;
    elseif img.wfilt==1
        img.res(i).sb(3).qstep = 4/wnorm_97(4,img.nwdec+2-i);
    end
    img.res(i).sb(3).sign = Ds;
    img.res(i).sb(3).coeff = uint16(floor(Dm/img.res(i).sb(3).qstep));
end

%%
% Now that wavelet coefficients are quantized and mapped onto a fixed
% number of bits, we can truly observe the scalability offered by the
% bit-plane progression used in JPEG 2000. To illustrate this, we choose a
% subband, let's say the HL subband of last resolution, and display its k most
% significant bit-planes for k=1,..,K where K=number of significant
% bit-planes for this subband.

sbcoeff = img.res(end).sb(1).coeff;
sbsign = img.res(end).sb(1).sign;
K = ceil(log2(max(double(sbcoeff(:)+1))));
h=figure;
for k=1:K
    subplot(2,ceil(K/2),k);
    nbp_discard = K-k;
    if nbp_discard
        mask_AND = bitshift(uint16(65535),nbp_discard); % 65535 = maximum dynamic with UINT16
        mask_OR = bitshift(uint16(1),nbp_discard-1);
        if nbp_discard==1
            mask_OR=0;
        end
        m_trunc = bitor(bitand(sbcoeff,mask_AND),mask_OR);
        m_trunc(sbcoeff<=pow2(nbp_discard)-1)=0;
        c_trunc = sbsign .* double(m_trunc);
    else
        c_trunc = sbsign .* double(sbcoeff);
    end
    imshow(c_trunc,[-max(abs(c_trunc(:))) max(abs(c_trunc(:)))]);
    title(sprintf('%d most sign. bp',k),'FontSize', 12);
end
suptitle(sprintf('HL-subband of the last resolution of image "%s" : \nprogressive bit-plane refinement',img.path));

if save_fig==1
    saveas(h,sprintf('figures/fig%d.fig',nfig));
    saveas(h,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end


%%
% As we see in the figure, wavelet coefficients are progressively refined,
% as the bit-planes (from the most to the least significant one) are
% included in the coefficient estimation. This principle will be used in
% the R-D allocation process, as explained in Section 3.


%% 2.2.2 Context Distribution Table (CDT) computation
% In this section, we compute the probability distribution of getting a '1'
% or a '0'. This distribution will be used when the entropy coding step is
% performed. Note that in JPEG 2000, this distribution is computed and
% adapted dynamically and independently for each code-block. In this
% experiment, for the sake of simplicity, we compute it "offline", directly
% on the whole image.

%%
% As explained above, we will compare two kinds of distribution. The first
% one is based on the frequency of occurence of one and zero symbols over
% the whole image. In the second case, the probabilities of binary random
% variables are estimated conditionally to their context. The relevance of
% such approach has been intuitively justified in Section XXXX.

%%
% Practically, the context of a bit corresponds to a set of state variables
% related to (i) the coefficient to whom the bit belongs, and (ii) its
% neighbouring coefficients. In our case, two state variables were used.
% (a) The "significant" status of a coefficient. Let us remind that a
% coefficient is said to become significant in a given bit-plane if a '1'
% bit is encountered for the first time for this coefficient (all other
% more significant bits were '0's). (b) The "first refinement" status of a
% coefficient. Among already significant coefficients, we will distinguish
% those that became significant in the previous (more significant)
% bit-plane.

%%
% In this work, we have considered 12 different contexts. They are
% presented in the preamble of function GET_CONTEXT and in Table 1 and 2.
% This is actually a subset of the contexts used in the JPEG 2000 standard,
% which uses 19 different contexts. 
% The nine first ones are used for not yet significant coefficients :
% depending on the kind of subband and on the amount and location of
% already significant coefficients surrounding the one to be coded, a
% different context will be used. The subband type is indeed taken into
% account as it influences the most probable location of already
% significant coefficients. The three last contexts are used to code bits
% from already significant coefficients. Again, the choice between these
% three contexts is done based on the amount of already significant
% coefficients around the one to be coded. In addition, the
% "first-refinement" status is also taken into account. The sign of each
% coefficient is introduced "as is" in the codestream and is not
% entropy-coded (in JPEG 2000, several contexts are defined specifically
% for the sign coding).

%%
% Let us first initialize the number of contexts and the Context
% Distribution Table (CDT). This vector stores the probability of having a
% '1' for each context. Last element of the vector is used to store the
% global probability of getting a '1' on the whole image.

global nctxt;
nctxt = 12;
CDT = zeros(nctxt+1,1);

%%
% As explained in Section XXX_THEORY_J2K_XXX, before being entropy-coded, subbands are
% divided in small entities called code-blocks. The width and height of a
% code-block are defined in the parameters above (must be power of 2 in the
% JPEG 2000 standard). Each code-block will then be entropy-coded
% separately, starting from the most significant bit-plane to the least
% significant one. 

%%
% In the following embedded loops, each subband from each resolution level
% is processed and divided in such code-blocks. Then, each code-block is
% analyzed. First, the number of all-zero most significant bit-planes is
% computed (NBPS field of each code-block). These bit-planes will not
% actually be coded : the number of all-zero bit-planes will simply be
% stored in the header of the packet containing the code-block coded
% bit-stream. Then, each code-block is analyzed through the ANALYZE_CB
% function. This function will process each bit-plane starting from the
% first most significant non-zero bit-plane and return two variables : (1)
% a matrix ctxt that gives two values for each context and each bit-plane :
% (i) the total number of bits and (ii) the number of '1' bits, (2) a
% vector DISTO that computes, for each bit-plane BP, the square error
% between the original code-block and the code-block whose least
% significant bit-planes are truncated, starting from and including
% bit-plane BP.
%
% Then, the DISTO values are adapted to reflect the square error in the
% pixel domain. To do so, they are multiplied by the squared value of the
% quantization stepsize and the squared norm of the corresponding synthesis
% filter.

CDT_tmp = zeros(2,nctxt);
for resno=1:numel(img.res)
    for sbno=1:numel(img.res(resno).sb)
        coeff = img.res(resno).sb(sbno).coeff;
        Ssign = img.res(resno).sb(sbno).sign;
        cbno=1;
        for y0=1:img.cbh:size(coeff,1)
            for x0=1:img.cbw:size(coeff,2)
                x1 = min(x0+img.cbw-1,size(coeff,2));
                y1 = min(y0+img.cbh-1,size(coeff,1));
                cb.coeff = coeff(y0:y1,x0:x1);
                cb.sign = Ssign(y0:y1,x0:x1);
                cb.nbps = ceil(log2(max(double(cb.coeff(:)+1))));
                if resno==1
                    cb.sbtype=sbno;
                else
                    cb.sbtype=sbno+1;
                end
                [cb.ctxt,cb.disto] = analyze_cb(cb);
                if img.wfilt==0
                    cb.disto = cb.disto.*((img.res(resno).sb(sbno).qstep)^2*...
                    (wnorm_53(cb.sbtype,img.nwdec+2-resno))^2);
                elseif img.wfilt==1
                    cb.disto = cb.disto.*((img.res(resno).sb(sbno).qstep)^2*...
                    (wnorm_97(cb.sbtype,img.nwdec+2-resno))^2);
                end
                
                img.res(resno).sb(sbno).cb(cbno)=cb;
                
                CDT_tmp=CDT_tmp+permute(sum(cb.ctxt,1),[3 2 1]);
                cbno=cbno+1;
            end
        end
    end
end

%%
% Once matrix CTXT has been computed for each code-block, we can easily
% determine the Context Distribution Table (CDT) that will be used for entropy
% coding. 

CDT(1:end-1) = (CDT_tmp(2,:)./CDT_tmp(1,:))';
CDT(end) = sum(CDT_tmp(2,:))/sum(CDT_tmp(1,:));

%%
% The CDT stores for each of the 12 contexts the probability to get a '1'
% bit according to the processed image. A 13th value stores the global
% probability to get a '1' bit, independently from the neighborhood of the
% coefficient. These values are presented in the following
% bar-graph.

h=figure;
title(sprintf('Context Distribution Table for image "%s"',img.path));
bar(CDT);
xlabel('Context number');
ylabel('P(1)');
axis([0 14 0 1]);
set(gca,'XTickLabel',{'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'Global P'});

if save_fig==1
    saveas(h,sprintf('figures/fig%d.fig',nfig));
    saveas(h,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end

%%
% The CDT is then saved so that it can be re-used when encoding other
% images.

[rext,rbase] = strtok(img.path(end:-1:1),'.');
CDTpath = [rbase(end:-1:2) '_CDT.mat'];
save(CDTpath,'CDT');

%% 2.2.3 R-D pairs computation based on CDT
% Once the CDT has been computed, we can use it to entropy-encode each
% code-block. In JPEG 2000, the entropy coder used is a binary arithmetic
% coder (MQ-coder, see~\cite{QCoder,normej2k}). In this experiment, we will
% limit ourselves to the computation of the rate that would be obtained if
% each code-block had been entropy-coded using this CDT. Indeed, the goal
% is not here to actually compress the image, but rather to assess the
% efficiency of context-based modeling and to compute the required data for
% the rate-allocation step.

%%
% Practically, for each subband, we progressively fill in a matrix RD :
% IMG.RES(RESNO).SB(SBNO).RD(CBNO,BPNO,i) where i=1 for the rate values and
% i=2 for the distortion values. Each pair of values (i=1,2) in the matrix
% RD gives therefore the amount of bits that will be needed to encode a
% given bit-plane from a given code-block in the processed subband (R), and
% the distortion reduction that this bit-plane bring when decoding the
% image (D). As all code-blocks do not have the same number of significant
% bit-planes in a subband, dimension 2 of matrix RD (counting the
% bit-planes) is taken equal to the maximum number of significant
% bit-planes among the code-blocks from the subband.
% 
% Note : 
%  - RESNO=1 for the smallest resolution.
%  - BPNO=1 corresponds to the most significant bit-plane.

%%
% If img.useownCDT = 0, we do not use the CDT computed on the processed
% image but the one specified in img.CDT2use. It should be noted that using
% a context distribution computed on an image whose content is very
% different from the one of the processed image does not change the
% performances drastically. This can be explained by the fact that the
% distribution is computed on all bit-planes of all subbands and that
% globally the distribution of each context does not depend that much on
% the content. In JPEG 2000 however, the distribution is re-initialized for
% each code-block and is dynamically adapted while encoding the bit-planes.
% In this case, the entropy coder uses a more adapted distribution and is
% therefore more efficient. As the goal here is more to highlight the
% improvement that can be drawn from the use of contexts compared to a non
% context-based solution, we keep a simple global static context
% distribution.

if img.useownCDT==0
    CDT=load(img.CDT2use);
    names = fieldnames(CDT);
    CDT=CDT.(names{1});
end


%%
% A matrix R_RES is also filled in. This matrix will contain, for each
% resolution level, 3 values : (1) The total number of bits that have to be
% compressed in this resolution level. (2) The total number of bits
% obtained when the entropy-coding step is done without contexts (global
% distribution on the image) (3) The total number of bits obtained when the
% entropy-coding step is done with contexts.
% 
% This matrix will then be used to compare the coding efficiency of the
% two different distributions (with and without contexts).

R_res = zeros(numel(img.res),3);

%%
% In these four embedded loops, we scan each bit-plane from each code-block
% in every subband from each resolution level. For each of them, we compute
% the rate RC that would be obtained after a context-based entropy coding,
% the rate RNC that would be obtained after a non-context-based entropy
% coding (i.e. with a global probability distribution on the whole image),
% and the distortion reduction brought by the decoding of the processed
% bit-plane (i.e. the difference between square error obtained with
% successive bit-planes).
%
% In particular, two functions are used inside the loops : (1) GET_RATE
% that computes the entropy of each bit-plane based on the given CDT. (2)
% GET_DISTO that computes the difference between square errors obtained with
% successive bit-planes.
%
% The matrix RD of each subband will then be used in the next section to
% find for each code-block the best truncation point, that is, the one that
% will minimize the distortion for a given global bit budget.

for resno=1:numel(img.res)
    for sbno=1:numel(img.res(resno).sb)
        sb = img.res(resno).sb(sbno);
        nbps = max([sb.cb(:).nbps]);
        ncb = numel(sb.cb);
        RD = zeros(ncb,nbps,2);
        for cbno=1:ncb
            ctxt=sb.cb(cbno).ctxt;
            coeff=sb.cb(cbno).coeff;
            disto=sb.cb(cbno).disto;
            offset=nbps-sb.cb(cbno).nbps;
            for bpno=1:sb.cb(cbno).nbps
                [rc rnc]=get_rate(ctxt,coeff,bpno,CDT);
                RD(cbno,offset+bpno,1)=rc;
                RD(cbno,offset+bpno,2)=get_disto(disto,bpno);
                R_res(resno,2)=R_res(resno,2)+rnc;
            end
            R_res(resno,1)=R_res(resno,1)+numel(sb.cb(cbno).coeff)*sb.cb(cbno).nbps...
                +numel(sb.cb(cbno).coeff(sb.cb(cbno).coeff>0));     
        end
        img.res(resno).sb(sbno).nbps = nbps;
        img.res(resno).sb(sbno).RD = RD;
        R_res(resno,3)=R_res(resno,3)+sum(sum(RD(:,:,1),1),2);
    end
end

%%
% Finally, we plot two figures comparing both approaches : with and without
% contexts. The first figure shows the compression ratio obtained for each
% resolution level. The total number of uncompressed bits taken to compute
% these values is the number of bits truly processed by the entropy coder,
% i.e. wavelet coefficients excluding the non-significant bit-planes. The
% second figure computes the global compression ratio, comparing the
% original number of bits in the pixel domain, with the rate obtained after
% the entropy coding step. As expected, the context-based approach is more
% efficient than the other one, as it more accurately estimates the
% probability of getting a '1' or a '0' depending on the coefficient
% location. In particular, in the graph plotting the compression ratio per
% resolution level, we see that the non-context-based entropy coding
% actually expands the LL-subband rather than compresses it. This is
% because the global probability distribution used in this case is very
% different than the one really observed for this resolution level. On the
% contrary, when using contexts, and even if those contexts do not take
% into account the resolution level, the entropy coder still achieves
% compression on the low frequencies.

g_c_eff=figure;
set(g_c_eff,'DefaultAxesFontSize',14);
bar(R_res(:,2:3)./repmat(R_res(:,1),1,2));
title('Lossless compression efficiency in the wavelet domain','FontSize', 14);
xlabel('Resolution number (1=LL, ...)','FontSize', 14);
ylabel('Compression ratio','FontSize', 14);
legend('No contexts','Context-based entropy coding');
axis([0 numel(img.res)+1 0 1.2]);
set(gca,'FontSize',12);

if save_fig==1
    saveas(g_c_eff,sprintf('figures/fig%d.fig',nfig));
    saveas(g_c_eff,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end

g_c_eff2=figure;
set(g_c_eff2,'DefaultAxesFontSize',14);
cratio(1,1) = sum(R_res(:,2))./(img.h*img.w*img.bdepth);
cratio(1,2) = sum(R_res(:,3))./(img.h*img.w*img.bdepth);
cratio(2,:) = 0;
bar(cratio);
title('Global lossless compression efficiency','FontSize', 14);
ylabel('Compression ratio','FontSize', 14);
legend('No contexts','Context-based entropy coding');
set(gca,'FontSize',12);
set(gca,'XTick',2);
axis([0.5 1.5 0 1]);

if save_fig==1
    saveas(g_c_eff2,sprintf('figures/fig%d.fig',nfig));
    saveas(g_c_eff2,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end


%% 2.3 Rate-Distortion optimal allocation
%
% In this section, we demonstrate the benefit of rate-distortion optimal bit allocation 
% when a target bit-budget has to be distributed across image blocks through the selection 
% of a quantization parameter within a discrete set of parameters corresponding to a discrete
% set of RD trade-offs. 
% Specifically, we consider an image whose wavelet coefficients have been split into a set 
% of codeblocks. The coefficients of each codeblock are encoded bitplane by bitplane, which 
% defines M distortion/bitbudget pairs (d_{ij},b_{ij}), corresponding to the approximation of the 
% coefficients of the $i^{th}$ codeblock by its $j^{th}$ most significant bitplanes, with 0 < j < M. 
%
% To illustrate the benefit of RD optimal bit allocation, we compare it to a naive allocation strategy 
% that simply assigns a constant number of bitplanes to all image codeblocks. In contrast, the RD optimal 
% strategy first computes the convex-hull RD points and select them in decreasing order of distortion 
% reduction per cost unit (see Section???). 
%

%% 2.3.1 Convex-hull approximation for each codeblock
%
% We assume the existence of a structure img.res(resno).sb(sbno).RD(i,j,k) that conveys, 
% for all image resolution and image subband, the incremental reduction of distortion and the 
% increase of bits corresponding to the addition of the $j^{th}$ bitplane to the definition of 
% codeblock $i$. Index $k$ is used to differenciate the cost in bits (k=1) from the decrease 
% in distortion (k=2).
%
% To implement the RD optimal bit allocation method, a 'hull structure' $HS$ is defined to collect 
% the convex-hull RD points of every codeblock in the image. As explained in Section????, for each 
% codeblock, those points are defined to be the subset of RD operating points that lie on the 
% lower convex-hull in a RD graph.  
%
% For each convex-hull RD point, the $HS$ structrure records the associated decrease of distortion 
% per bit unit, and a set of parameters characterizing the operating point. 
% Formally, each row of the table collecting convex-hull RD points has the form
% [gain, res_index, sb_index, cdblock_index,nbplanes, deltaR, deltaD], in which 'gain'
% refer to the reduction of distortion per unit of rate measured in
% comparison to the previous convex-hull RD point for the codeblock, while
% other parameters respectively define the codeblock at hand, the number
% of bitplanes associated to the operating point, and the increment (decrement) 
% of rate (distortion) compared to the previous convex-hull RD point of the codeblock.


nbHS = 0; % Number of rows of the table collecting all convex-hull RD operating points.

% Dealing with low-pass component.
for k=1:size(img.res(1).sb(1).RD,1)
            RD = img.res(1).sb(1).RD(k,:,:);
            lhbp=0;                          % Number of bitplanes corresponding to the last operating point found on the convex-hull.
            gain_max=1;                      % Just to make sure we enter in the while loop.
            
            % Add convex-hull RD points in increasing order of bitbudget,
            % until the operating point with minimal distortion has been reached.
            while (gain_max>0 && lhbp<size(RD,2))
            nhbp=lhbp;
            gain_max=0;
            
            % Starting from the last convex-hull RD point, we go through all subsequent (with larger rate) 
            % operating point of the block, searching for the one with maximal distortion reduction per bit. 
            % This operating point defines the next convex-hull RD point.
            for bp=(lhbp+1):size(RD,2)      
                if sum( img.res(1).sb(1).RD(k,(lhbp+1):bp,1) ) >0
                    gain=sum( img.res(1).sb(1).RD(k,(lhbp+1):bp,2) ) / sum( img.res(1).sb(1).RD(k,(lhbp+1):bp,1) );
                else gain = 0;
                end
                if gain > gain_max
                gain_max=gain;
                nhbp=bp;
                end
            end
            nbHS=nbHS+1;
            deltaR = sum( img.res(1).sb(1).RD(k,(lhbp+1):nhbp,1) );
            deltaD = sum( img.res(1).sb(1).RD(k,(lhbp+1):nhbp,2) );
            HS(nbHS,:)=[gain_max,1,1,k,nhbp, deltaR, deltaD];
            lhbp=nhbp;
            end
end

%%
% Same code as above, but dealing with all details subbands 
% instead of LL component.

% Dealing with all resolutions.
for i=2:img.nwdec+1
    % Dealing with all subbands.
    for j=1:3
        for k=1:size(img.res(i).sb(j).RD,1)
            RD = img.res(i).sb(j).RD(k,:,:);
            lhbp=0;                          
            gain_max=1;                                 
            while (gain_max>0 && lhbp<size(RD,2))
            nhbp=lhbp;
            gain_max=0;            
            for bp=(lhbp+1):size(RD,2)
                if sum( img.res(i).sb(j).RD(k,(lhbp+1):bp,1) )>0
                    gain=sum( img.res(i).sb(j).RD(k,(lhbp+1):bp,2) ) / sum( img.res(i).sb(j).RD(k,(lhbp+1):bp,1) );
                else gain = 0;
                end
                if gain > gain_max
                gain_max=gain;
                nhbp=bp;
                end
            end
            nbHS=nbHS+1;
            deltaR = sum( img.res(i).sb(j).RD(k,(lhbp+1):nhbp,1) );
            deltaD = sum( img.res(i).sb(j).RD(k,(lhbp+1):nhbp,2) );
            HS(nbHS,:)=[gain_max,i,j,k,nhbp, deltaR, deltaD];
            lhbp=nhbp;
            end
        end
    end
end



%% 2.3.2 RD curves computation
%
% Now that the hull structure has been computed, we can compare the naive 
% and RD optimal allocation strategies.
% The naive startegy assigns a constant number of bitplanes to each codeblock.
% For comparison purposes, the RD optimal allocation targets the same 
% bit-budget as the one achieved with naive allocation. 


%% Naive allocation
%
% Constant number of bitplanes assigned by naive allocation 
% is chosen from 1 to NBP.
NBP = 7;    

for n=1:NBP
% The global image (R,D) point corresponding to a constant number n of 
% bitplanes allocated to each codeblock is computed.
% NR and ND variables record the image rate/distortion trade-offs obtained with
% naive allocation.
    NR(n)=0;
    ND(n)=0;
    tmp_n = min(n, size(img.res(1).sb(1).RD, 2));
    NR(n) = NR(n) + sum( sum(img.res(1).sb(1).RD(:,1:tmp_n,1), 2), 1) ;
    ND(n) = ND(n) + sum( sum(img.res(1).sb(1).RD(:,:,2), 2), 1) - sum( sum(img.res(1).sb(1).RD(:,1:tmp_n,2), 2), 1);
    for i=2:img.nwdec+1
        for j=1:3
            tmp_n = min(n, size(img.res(i).sb(j).RD, 2));
            NR(n) = NR(n) + sum( sum(img.res(i).sb(j).RD(:,1:tmp_n,1), 2), 1) ;
            ND(n) = ND(n) + sum( sum(img.res(i).sb(j).RD(:,:,2), 2), 1) - sum( sum(img.res(i).sb(j).RD(:,1:tmp_n,2), 2), 1) ;        
        end
    end
end

NPSNR = -10*log10(ND/(255*255*img.h*img.w));

%% RD optimal allocation
%
% Once the HS structure has been defined, the RD optimal allocation strategy simply consists 
% in selecting the convex-hull RD points in decreasing order of distortion reduction per 
% cost unit, independently of the codeblock index. 
%
% Hence, the first step consists in sorting the HS structure.

ascendingHS = sortrows(HS);

%%
% We can then easily compute the RD optimal allocation for any target 
% bitbudget by reading the rows of 'ascendingHS' in decreasing order 
% of index, until the targettted bitbudget has been reached. 
%
% While computing a RD optimal allocation, we also want to keep track of the bitplanes 
% assigned to each codeblock so as to reconstruct the corresponding image at any time. 
% For this reason, the number of bitplanes assigned to the codeblocks are stored in the
% img.res(i).sb(j).nbrplanes array, initialized as follows:

img.res(1).sb(1).nbrplanes = zeros( 1, size(img.res(1).sb(1).RD,1) ); 
for i=2:img.nwdec+1
   img.res(i).sb(1).nbrplanes = zeros( 1, size(img.res(i).sb(1).RD,1) );
   img.res(i).sb(2).nbrplanes = zeros( 1, size(img.res(i).sb(2).RD,1) );
   img.res(i).sb(3).nbrplanes = zeros( 1, size(img.res(i).sb(3).RD,1) );
end

%%
% The images that are reconstructed to visually compare the naive and RD optimal
% allocation strategies are defined based on the NB parameter. 
% Specifically, the two reconstructed images are the ones obtained with naive 
% and RD optimal allocation for the same rate, defined to be equal to the one 
% obtained when NB bitplanes are assigned to each codeblock. 

NB = 2;     % The value of this parameter has been selected arbitrarily.


%%
% Rtmp and Dtmp are temporary variable that accumulate the rate and 
% decrement the image distortion along the allocation process.

Rtmp=0;   

%%
% Dtmp is initialized to the total image distortion measured when no
% bitplane is decoded.

Dtmp = sum( sum(img.res(1).sb(1).RD(:,:,2), 2), 1);
for i=2:img.nwdec+1
    for j=1:3
            Dtmp = Dtmp + sum( sum(img.res(i).sb(j).RD(:,:,2), 2), 1) ;        
    end
end

%%
% The RD optimal allocation process simply selects rows of 'ascendingHS' 
% in decreasing order of index.

indHS = nbHS;       % Index of the next convex-hull RD points to consider. 

%%
% We consider the global image RD points with rate equal to the rate 
% obtained with naive allocation of n bitplanes to each codeblock 
% (1<= n <= NBP).
% In the code below, the OR(n) and OD(n) variables are defined to record 
% the image rate distortion trade-offs obtained with RD optimal allocation, 
% for target rates corresponding to the rates achieved when assigning n 
% bitplanes to all image codeblock.

for n=1:NBP         
    while (Rtmp < NR(n) && indHS>0)
       Rtmp = Rtmp + ascendingHS(indHS,6);
       Dtmp = Dtmp - ascendingHS(indHS,7);
       
        % The purpose of this part of the code is to record the number of bitplanes decoded 
        % for each codeblock until the rate goes beyond the rate achieved when allocating 
        % NB bitplanes to each codeblock.                  
       if(n<=NB)      
          img.res(ascendingHS(indHS,2)).sb(ascendingHS(indHS,3)).nbrplanes(ascendingHS(indHS,4)) = ascendingHS(indHS,5); 
       end
       
       indHS = indHS - 1;
    end
    OR(n) = Rtmp;
    OD(n) = Dtmp;
end

%%
% Translation of distortion measures in PSNR values.

OPSNR = -10*log10(OD/(255*255*img.h*img.w));    


%% 2.3.3 Plotting and comparing the allocation methods.
%
% Here, we just plot PSNR vs Bit Budget graphs, and compare RD optimal and naive
% allocation methods.

g=figure;
set(g,'DefaultAxesFontSize',14);
plot(NR/1000,NPSNR,'b-.o',OR/1000,OPSNR,'r-^','LineWidth',1, 'MarkerSize', 4);
title('Rate-Distortion curves');
xlabel('Bit budget (kbits)','FontSize', 14);
ylabel('PSNR (dB)','FontSize', 14);
legend('Naive', 'RD Optimal');
grid on;
axis([0 1200 24 64]);

if save_fig==1
    saveas(g,sprintf('figures/fig%d.fig',nfig));
    saveas(g,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end


%% 2.4 Image reconstruction and display
%
% Here, Matlab code is provided to reconstruct the image based on the number 
% of bitplanes allocated to each codeblock. It permits to compare both 
% allocation methods from a perceptual point of view.  

%% 2.4.1 Code-blocks truncation based on RD results & inverse quantization
% We use the NBRPLANES field from each subband to truncate the
% code-blocks according to the optimal allocation computed in the previous
% Section. As truncation is equivalent to a quantization process, we drop
% the least significant bit-planes and change the coefficient value so that
% it lieas at the center of the quantization stepsize. This is done by
% setting to one the bit immediately following the last bit that has been
% kept.
%
% Then, an inverse quantization is applied on the truncated coefficients.
% This is simply done by multiplying the coefficients by the quantization
% step size used during the encoding.

NC_R = zeros(1,img.h*img.w);    % Reconstructed coefficients in case of naive allocation.
OC_R = zeros(1,img.h*img.w);    % Reconstructed coefficients in case of optimal allocation.
C_R_offset = 1;
for resno=1:numel(img.res)
    for sbno=1:numel(img.res(resno).sb)
        sb = img.res(resno).sb(sbno);
        ncb = numel(sb.cb);
        [h,w] = size(sb.coeff);
        NcoeffR = uint16(zeros(h,w));
        OcoeffR = uint16(zeros(h,w));
        h_offset = 1;
        w_offset = 1;
        for cbno=1:ncb
            bp_offset=sb.nbps-sb.cb(cbno).nbps;
            coeff=sb.cb(cbno).coeff;
            Nnbp=min(NB,size(sb.RD,2)) -bp_offset;              
            Onbp=sb.nbrplanes(cbno) -bp_offset;
            
            if Nnbp <= 0
                NcoeffR_tmp = uint16(zeros(size(coeff)));
            else
                nbp_discard = sb.cb(cbno).nbps-Nnbp;
                mask_AND = bitshift(uint16(65535),nbp_discard); % 65535 = maximum dynamic with UINT16
                mask_OR = bitshift(uint16(1),nbp_discard-1);
                if nbp_discard==1
                    mask_OR=0;
                end
                NcoeffR_tmp = bitor(bitand(coeff,mask_AND),mask_OR);
                NcoeffR_tmp(coeff<=pow2(nbp_discard)-1)=0;
            end
            
            if Onbp <= 0
                OcoeffR_tmp = uint16(zeros(size(coeff)));
            else
                nbp_discard = sb.cb(cbno).nbps-Onbp;
                mask_AND = bitshift(uint16(65535),nbp_discard); % 65535 = maximum dynamic with UINT16
                mask_OR = bitshift(uint16(1),nbp_discard-1);
                if nbp_discard==1
                    mask_OR=0;
                end
                OcoeffR_tmp = bitor(bitand(coeff,mask_AND),mask_OR);
                OcoeffR_tmp(coeff<=pow2(nbp_discard)-1)=0;
            end

            
            NcoeffR(h_offset:h_offset+size(coeff,1)-1,w_offset:w_offset+size(coeff,2)-1)=NcoeffR_tmp;
            OcoeffR(h_offset:h_offset+size(coeff,1)-1,w_offset:w_offset+size(coeff,2)-1)=OcoeffR_tmp;
            
            if (w_offset+size(coeff,2)-1)==w
                w_offset=1;
                h_offset=h_offset+size(coeff,1);
            else
                w_offset=w_offset+size(coeff,2);
            end
        end
        NcoeffRQ = sb.sign.*(double(NcoeffR)*sb.qstep);
        OcoeffRQ = sb.sign.*(double(OcoeffR)*sb.qstep);
        NC_R(C_R_offset:C_R_offset+numel(NcoeffRQ)-1)=NcoeffRQ(:)';
        OC_R(C_R_offset:C_R_offset+numel(OcoeffRQ)-1)=OcoeffRQ(:)';
        C_R_offset = C_R_offset+numel(NcoeffRQ);
    end
end

%% 2.4.2 Inverse Discrete Wavelet Transform
% After truncation and inverse quantization, an inverse DWT is applied to
% get back in the pixel domain.

if img.wfilt==0
    NX_R = waverec2(NC_R,S,lo_53_R,hi_53_R);                  % dyadic decomposition
    OX_R = waverec2(OC_R,S,lo_53_R,hi_53_R);                  % dyadic decomposition
elseif img.wfilt==1
    NX_R = waverec2(NC_R,S,lo_97_R,hi_97_R);                  % dyadic decomposition
    OX_R = waverec2(OC_R,S,lo_97_R,hi_97_R);                  % dyadic decomposition
else
    error('wavelet filter not recognized');
end

%% 2.4.3 Post-processing
% An inverse DC level shifting is finally applied on the reconstructed
% pixel values.

NX_R = NX_R + pow2(img.bdepth-1); % inverse DC LEVEL SHIFTING
OX_R = OX_R + pow2(img.bdepth-1); % inverse DC LEVEL SHIFTING


%% 2.4.4 Reconstructed image display

h=figure;
imshow(NX_R,[0 255]);
title(sprintf('Reconstructed image of %s at %4.2f kbits, Naive allocation.',img.path, NR(NB)/1000 ));

if save_fig==1
    saveas(h,sprintf('figures/fig%d.fig',nfig));
    saveas(h,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end

h=figure;
imshow(OX_R,[0 255]);
title(sprintf('Reconstructed image of %s at %4.2f kbits, RD optimal allocation.',img.path, OR(NB)/1000 ));

if save_fig==1
    saveas(h,sprintf('figures/fig%d.fig',nfig));
    saveas(h,sprintf('figures/fig%d.pdf',nfig));
    nfig=nfig+1;
end






