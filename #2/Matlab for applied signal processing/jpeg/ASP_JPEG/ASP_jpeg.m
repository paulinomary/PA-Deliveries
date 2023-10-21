%% Chapter 8 - How are digital images compressed in the web?
%
% F. Marques, M. Menezes, J. Ruiz 2008
%
% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%
% This file uses the image toolbox of MATLAB.

%% Proof of concepts
% In this Section we are going to illustrate the main concepts behind the
% JPEG standard. These concepts are applied as well in video compression in
% the so-called hybrid approach which is the basis for standards such as
% MPEG1 or MPEG2 and, therefore, these concepts will be revisited in
% Chapter 9. Note that, in these two Chapters, we will not concentrate on
% the specificities of any of these standards but on the common, generic
% tools and concepts that they use. Moreover, in order to illustrate these
% concepts, we will work mainly with grey level images and explain the
% extensions that should be introduced in the proposed systems to work with
% color images. Let us start by reading a grey level image. In our case,
% this is the first image of a video sequence that, in this Chapter, we
% will use it as an example to study the case of still (isolated) image
% coding and of the exploitation of the spatial redundancy. Moreover, in
% Chapter 9, it will be used as the initial point for the study of video
% sequence coding and of the exploitation of both spatial and temporal
% redundancies. First we will study how the spatial redundancy can be
% exploited to reduce the number of bits needed to represent the image. Let
% us load file 'table_000_g.bmp' from disk and show its pixel values (Fig.
% 8.11). This file contains a gray image corresponding to frame #000 of the
% Table Tennis sequence with size 176x144 pixels  stored in non compressed
% form.

iptsetpref('ImshowBorder', 'tight');format compact;

table = imread('seqs/table/table_000_g.bmp');
imshow(table);

%% Block Image Transformation

%% Image block processing

% Spatial redundancy is exploited by performing a local analysis of the
% image; that is, by dividing the image into non-overlapping square blocks.
% Then, the information contained in each of these blocks will be processed
% separately. Typically, images are partitioned into blocks of 8x8 pixels
% (Fig. 8.12).

addgridtofigure(size(table),[8 8]);
title('block processing for image Table tennis');

%% Mirror padding

% Given that generic images (images of any size) may not allow an exact
% division in 8x8 pixel blocks, information in the right and low boundary
% blocks is usually padded before partitioning the image. A common way to
% pad this information is by mirroring it (which reduces the transitions
% introduced when padding with, for instance, zero values). Let us see this
% case with a new version of the previous Table Tennis image where 4 rows
% and 4 columns have been removed from the bottom and right sides of the
% original image (Fig. 8.13).

table_crop = table(1:end-4,1:end-4);
table_padded = padarray(table_crop,[4,4],'symmetric','post');
imshow(table_padded);
title('Mirror padded table\_crop image');

% If we analyze one of the mirrored blocks, the effect of the mirroring can
% be easily seen (Fig. 8.14). For instance, the block in the first row and
% last column shows how the mirror padding is created along the vertical
% axis. Note as well how this padding preserves the original block
% contours.

imshow(kron(table_padded(1:8,end-7:end),uint8(ones(16))));
title('Mirror padded block');

%% DCT block transformation

% Now, we are going to analyze how our basic coding unit (that is, a
% generic 8x8 pixel block which, in the sequel, we will refer to as a
% block) is processed. Among all the different transforms that can be
% applied on an image block to obtain a set of less correlated coefficients
% in the transformed domain, the Discrete Cosine Transform (DCT) has been
% shown to present very good properties (see the discussion in the previous
% Section).

% Let us analyze one block from the Table Tennis image to illustrate the
% whole transform process. First we select a given block from the image,
% for example, the block situated 12 blocks from the top and 13 from the
% left:

posr = 12; posc = 13;
f = table*0.2;
f(8*(posr-1)+1:8*(posr-1)+8,8*(posc-1)+1:8*(posc-1)+8) = ...
    table(8*(posr-1)+1:8*(posr-1)+8,8*(posc-1)+1:8*(posc-1)+8);
imshow(f);
%addgridtofigure(size(table),[8 8]);

% The following image shows a magnification of the selected block. Also,
% the corresponding matrix of pixel values is shown:

b = table(8*(posr-1)+1:8*(posr-1)+8,8*(posc-1)+1:8*(posc-1)+8); 
imshow(kron(b,uint8(ones(16))));
b,

% The MATLAB function dct2 performs the DCT of a given image block. The
% transformed matrix (that is, the matrix containing the DCT coefficient
% values) for the block under study is

d = dct2(double(b)-128),

% where as explained in the previous Section a value of 128 has been
% subtracted to all pixels of block b to create a signal with zero mean and
% obtain a lower energy DC coefficient. We can observe the capability of
% the DCT to compact the energy of the signal with respect to the energy
% distribution in the original block. Fig. 8.17 shows the absolute value of
% the DCT coefficients (sorted in descending order):

v = sort(abs(d(:)),1,'descend');
figure,plot(v),grid on; 
xlabel('Coefficient number');
ylabel('Absolute value');
title('sorted DCT coefficients');

% If we plot the percentage of accumulated energy for these coefficients
% (see Fig. 8.18), it can be noted that, for this specific block, 95% of
% the energy of the signal in the transform domain is contained in the 13
% coefficients with highest energy.

dct_energy_distribution = cumsum(v.*v)/sum(v.*v);
plot(dct_energy_distribution),grid on
xlabel('Coefficient number'); 
ylabel('Accumulated % energy');
title('Accumulated percentage of energy for the sorted DCT coefficients');

% However, in the original signal (the block b), the same percentage of
% energy is reached when adding up the contributions of the 47 pixels with
% the highest energy (as shown if Fig. 8.19). 

bv = sort(double(b(:)),1,'descend');
ima_energy_distribution = cumsum(bv.*bv)/sum(bv.*bv);
plot(ima_energy_distribution),grid on
xlabel('Coefficient number'); 
ylabel('Accumulated % energy');
title('Accumulated percentage of energy for the sorted pixels');

% Moreover note that, roughly speaking, coefficients in the DCT domain with
% higher energy are gathered around the left-top corner of matrix d (lower
% frequency coefficients). This is a common behavior for all natural image
% block DCT coefficient matrices, which will be further exploited in the
% coding process when including the quantization step.

% The DCT transformation is invertible and we can recover the initial pixel
% values of the original block b using the inverse DCT transformation
% (function idct2 in MATLAB):

br = uint8(idct2(d)+128);
b - br,


%% Elimination of DCT coefficients

% Now, we can see the effect of zeroing some of the DCT coefficients in the
% reconstruction of the original block. In Fig. 8.20 we present the result
% of reconstructing the image block using the first N = (1, 4, 8, 16, 32,
% 64) coefficients in the zigzag scan. This strategy to select coefficients
% is often referred to as “zonal coding”.

v = zigzag8(d);
dct_energy_distribution = cumsum(v.*v)/sum(v.*v);
for N=[1 4 8 16 32 48],
    ve = [v(1:N),zeros(1,64-N)];
    dr = izigzag8(ve);
    br = uint8(idct2(dr)+128);
    figure,imshow(kron(br,uint8(ones(16))));
    disp(['N = ' int2str(N) ' (~' int2str(round(dct_energy_distribution(N)*100)) '%)']);
end;

% In turn, in Fig. 8.21 we present the result of using the same number of
% coefficients as in the example but, in this case, the N coefficients with
% highest energy are selected directly. This strategy to select
% coefficients is often referred to as “threshold coding”. Below each
% image, the percentage of energy corresponding to the N selected
% coefficients is presented.

d2 = d(:);
[v,i] = sort(abs(d2),1,'descend');
v = d2(i); ii(i) = 1:64;
dct_energy_distribution = cumsum(v.*v)/sum(v.*v);
for N=[1 4 8 16 32 48],
    ve = [v(1:N);zeros(64-N,1)];
    veu = ve(ii);
    dr2 = reshape(veu,8,8);
    br = uint8(idct2(dr2)+128);
    imshow(kron(br,uint8(ones(16))));
    disp(['N = ' int2str(N) ' (~' int2str(round(dct_energy_distribution(N)*100)) '%)']);
end;

%% Complete image block coding

% In order to draw more general conclusions, we are going to analyze some
% of these features in a whole image. Matrix Ed presents the energy
% distribution per DCT coefficient averaged among all the blocks in the
% image. As it can be seen, the previous behavior is preserved; that is,
% coefficients with largest energy are those of lowest frequency (top-left
% matrix corner).

D = blkproc(double(table)-128,[8 8],@dct2);
Ed = blkmean(D.*D,[8 8]);
Ed = round(Ed/sum(sum(Ed))*100),

% The effect of dropping some DCT coefficients when reconstructing the
% complete image (the so-called Block Effect) can be observed in the
% following two Figures. Since an image independent block partition has
% been imposed and blocks have been separately coded, block boundaries are
% noticeable in the decoded image. This effect is more remarkable when
% fewer coefficients are used. In Fig. 8.22, we present six different
% versions of the original Table Tennis image: keeping N = (1, 4, 8, 16,
% 32, 64) coefficients at each block, respectively.

n = 1;
dct_total_energy = sum(sum(D.^2));
for N=[1 4 8 16 32 64],
    Dk = blkproc(D,[8 8],@coeffs_keep_zigzag,N);
    dct_energy = sum(sum(Dk.^2));
    tabler = uint8(blkproc(Dk,[8 8],@idct2)+128);
    figure,imshow(tabler);
    disp([char(96+n) ') N = ' int2str(N) ' (~' int2str(round(dct_energy/dct_total_energy*100)) '%)']);
    n = n + 1;
end;

% In Fig. 8.23 we present the result of using the same number of
% coefficients as in the previous example, but, in this case, those
% coefficients with higher energy are selected (“threshold coding”). For
% each number of coefficients N, the % of energy selected is shown.

n = 1;
dct_total_energy = sum(sum(D.^2));
for N=[1 4 8 16 32 64],
    Dk = blkproc(D,[8 8],@coeffs_keep_higher,N);
    dct_energy = sum(sum(Dk.^2));
    tabler = uint8(blkproc(Dk,[8 8],@idct2)+128);
    figure,imshow(tabler);
    disp([char(96+n) ') N = ' int2str(N) ' (~' int2str(round(dct_energy/dct_total_energy*100)) '%)']);
    n = n + 1;
end;

% The analysis of the previous examples leads to the following comments: 

% * Increasing the same amount of coefficients/energy in all blocks does
% not translate into the same increase of visual quality for all blocks.
% This is due to the fact that not all DCT coefficients are equally
% relevant for the human visual system. This concept will be further
% exploited in the sequel to adapt the quantization of the DCT coefficients
% to the final human observer.

% * In order to obtain a homogenous image quality (for instance, a given
% global PSNR), a different amount of coefficients/energy (and, eventually,
% bits) can be assigned to the various blocks. Therefore, given a bit
% budget for a specific image, those bits can be shared among the blocks in
% a content dependent manner. That is, more/fewer bits are used for those
% blocks containing more/less complex information. This concept is the
% basis for the Rate-Distortion Theory that will be further addressed in
% Chapter 10.



%% DCT quantization

% As commented in the previous Section, DCT coefficients need to be
% quantized in order to generate a discrete source whose symbols can be
% afterwards entropy encoded. In the current case, a separated scalar,
% uniform quantization is used for each DCT coefficient.

% The different quantization of each coefficient is represented by the
% so-called Quantization Table (Q Table). This is a matrix containing at
% position (i,j) the value to which the DCT coefficient (i,j) should be
% compared to generate its associated quantization index . A possible Q
% Table was proposed by Lohscheller  (Lohscheller 1984) and it is presented
% as matrix Q. As it can be seen, roughly speaking, the higher the DCT
% coefficient frequency, the larger the associated value in the Q Table.
% This table is  the results of psychovisual experiments that have
% determined the relative relevance of the different DCT coefficients 
% for the human visual system.

Q = jpegsteps,

% There is a difference in the way in which DCT coefficients are commonly
% indexed and their MATLAB implementation. DCT coefficients of an 8x8 pixel
% block are usually indexed from 0 to 7 in both dimensions with the first
% index indicating the desired column and the second the desired row. In
% turn, the coordinates of a MATLAB matrix are indexed from 1 to 8 in both
% dimensions with the first index indicating the desired row and the second
% the desired column.
  
% Lohscheller Tables have been adopted in the JPEG and MPEG standards as
% default proposals, although, in these standards, different tables can be
% set by the encoder and included in the bit stream. In the case of JPEG,
% there is a slight varia-tion in the Q(1,1) value since, due to
% implementation features, the implemented DCT transform is not unitary

  % The comparison between a given DCT coefficient and its associated value
% in the Q Table is carried out as expressed in the following code in
% which, for illustration purposes, we are using the same block we have
% used in the previous subsection (see Fig. 8.16): 

dq = round( d./Q ),

% For example, a coefficient d(3,1) smaller than 5 implies that the
% coefficient is zeroed after quantization (Q(3,1) = 10) whereas
% coefficient d(6,8) has to be smaller than 50 in order to be zeroed
% (Q(6,8) = 100).

% Different coding qualities can be obtained by multiplying the values in
% the Q Table by a given constant, k, producing a new Q Table: kQ. The
% following MATLAB code presents the quantization of the block b using k =
% 3 and k = 5.

k=3; 
k.*Q,round(d./k./Q),

k=5; 
k.*Q,round(d./k./Q),

% In turn, Fig. 8.24 presents the reconstruction of the Table Tennis image
% when quantizing the DCT coefficients of its blocks using k = 1, k = 3 and
% k = 5. In our implementation, the MATLAB functions quantizedct8 and
% iquantizedct8 perform the (direct and inverse) quantization of the DCT
% coefficients of a given matrix. Moreover, they are responsible of
% clipping the DC value between [0, 255] and AC values between [-128, 127]. 

n = 1;
tabledct = blkproc(table, [8 8], @dct2);
for k=[1 3 5],
    tabledctq = blkproc(tabledct, [8 8], @quantizedct8, k.*Q);
    tabledctqi = blkproc(tabledctq, [8 8], @iquantizedct8, k.*Q);
    tabler = uint8(blkproc(tabledctqi, [8 8], @idct2));
    figure,imshow(tabler);
    disp(['k = ' int2str(k)]);
    n = n + 1;
end;

% Note that, in the proposed implementation, the value Q(1,1) = 8 is used
% and there is no subtraction of 128 to the pixel values before computing
% the DCT trans-form. Therefore, the DC coefficient can be directly
% quantized in the range [0, 255].

%% Spatial decorrelation between blocks

% So far, we have decorrelated the information in the image only using the
% DCT and this has been done grouping the data into blocks of 8x8 pixels.
% Since this grouping is independent of the image content, it is likely
% that collocated transformed coefficients in neighbor blocks will still be
% correlated. If such a correlation is present, an additional
% transformation can be applied on the DCT coefficients to further reduce
% it. To illustrate this correlation, the next Figures present the
% 2-Dimensional distribution of collocated DCT coefficients in consecutive
% blocks. To ensure that the blocks that are compared are neighbors, blocks
% are ordered following the zigzag scan.

% In order to create more meaningful 2-Dimensional histograms (that is, to
% have more data), we use a 352x288 pixels image. This image format is
% usual known as Common Intermediate Format (CIF). Moreover, we have
% selected the first image of the Foreman sequence instead of the Table
% Tennis one for this specific example. The Table Tennis image has a large
% number of blocks with similar DC values (all blocks corresponding to the
% wall in the background), which prevents from correctly illustrating the
% spatial decorrelation effect. On the contrary, the Foreman sequence (see
% Fig. 8.25) presents more variation of the DC values as blocks are less
% homogeneous. Fig. 8.26, Fig. 8.27, and Fig. 8.28 show the distribution of
% the values of DCT coefficients in positions (1,1), (2,1) and (4,4) for
% consecutive blocks in the zigzag scan, respectively.

k = 1;
im = imread('seqs/foreman/foreman_cif_000_g.bmp');
imdct = blkproc(im,[8 8],@dct2);
imdctq = blkproc(imdct,[8 8],@quantizedct8,k*Q);
pos = {[1 1],[2 1],[1 3]};
for p = 1:length(pos),
    coeffsm = imdct(pos{p}(1):8:end,pos{p}(2):8:end);
    coeffs = zigzag(coeffsm);
    figure,hold on,grid
    for i=2:length(coeffs),
        plot(coeffs(i-1),coeffs(i),'o');
    end;
    if (p==1),axis([0 2000 0 2000]);end;
    if (p==2),axis([-500 500 -500 500]);end;
    if (p==3),axis([-200 200 -200 200]);end;
    xlabel('value DCT coeff i');
    ylabel('value DCT coeff i-1');
end;

% Note that the highest correlation appears in the case of the DC
% coefficient (at position (1,1)) since the points of the distribution are
% placed near the diagonal. This correlation is commonly exploited by using
% a prediction step to represent the DC values of consecutive blocks. This
% way, quantized DC coefficients are losslessly coded using a Differential
% Pulse Code Modulation (DPCM) technique. In the proposed implementation,
% the DC coefficient is predicted by the mean value of the DC coefficients
% of the 3 blocks above and the block on the left (Fig. 8.29).

coeffsdc = imdctq(1:8:end,1:8:end);
h = [-0.25 -0.25 -0.25; -0.25 1 0; 0 0 0];
coeffsdc_decorrelated = filter2(h, coeffsdc);
coeffs = zigzag(coeffsdc_decorrelated);
figure,hold on,grid
for i=2:length(coeffs),
    plot(coeffs(i-1),coeffs(i),'o');
end;
axis([-1000 1000 -1000 1000])
xlabel('value of coefficient i');
ylabel('value of coefficient i-1');


%% Entropy coding

% At this moment, we have two different sets of quantization indices or
% symbols: those related to the prediction error of the decorrelated
% quantized DC coefficients and those associated to the zigzag ordered AC
% coefficients.

% The JPEG standard specifies that these coefficients are entropy encoded
% through either Huffman or Arithmetic encoding. Commonly, Huffman encoding
% is used and, thus, this is the approach that we discuss in this Section. 

% Let us see the entropy coding in the proposed MATLAB implementation where
% Huffman coding is used. Let us start as an example with the quantized DCT
% coefficients of our selected block b (using k = 1).

k = 1;
dq = round(dct2(b)./k./Q),
 
% As explained earlier in this Chapter, the DC coefficient is decorrelated
% with previous DC values and encoded separately. In the case of the AC
% coefficients, as it can be seen in the example above, there are many AC
% coefficients with value 0 (even with k = 1). Therefore, it is important
% to precede the Huffman encoding of AC coefficients with a runlength
% encoder. Again, the AC coefficients are zig-zag scanned before the
% runlength encoder:

dc = dq(1,1);
ac = zigzag(dq); ac = ac(2:end);
[nz vv] = runlength(ac,63);

% The MATLAB function encodemessage is responsible for translating the
% quantized DC coefficients and the runlength AC coefficients into symbols
% using Huffman codes. The MATLAB file "codes/imcodes.mat" includes the
% previously generated Huffman codes (as seen in the previous Section)
% corresponding for each DC coefficient (in variable dcdcode) and each AC
% coefficient (variables nzcode and vvcode). In our example, the final list
% of symbols for the DC coefficient (in this case, as we analyze an
% isolated block, the DC coefficient has not been decorrelated) is:

load('codes/imcodes.mat');
codeddcd = encodemessage(dc+256', dcdcode),

% As the decorrelated DC values have a dynamic range between [-255, 255] a
% value of 256 is added to form the index to access Huffman codes stored in
% variable dcdcode.

% The rulength encoding of AC coefficients is (again, 2 is added to the
% runlength run values to create indexes as the dynamic range is [-1, 62]
% and 129 is added to runlength values with dynamic range [-128, 128]):

codednz = encodemessage(nz(:)' + 2, nzcode),
codedvv = encodemessage(vv(:)' + 129, vvcode),

% The following code implements the entropy coding for all blocks of the
% input image.

tabledct  = blkproc(table, [8 8], @dct2);
tabledctq = blkproc(tabledct, [8 8], @quantizedct8, k.*Q);

zzdctq = double(blkproc(tabledctq, [8 8], @zigzag8));

% Separate DC from AC components of the DCT.
[zznrows zzncolumns] = size(zzdctq);
mask = true([zznrows zzncolumns]);
mask(:, 1 : 64 : end) = false;
dc = zzdctq(~mask);
ac = reshape(zzdctq(mask), ...
    [zznrows, zzncolumns - zzncolumns / 64])' - 128;

% Decorrelate DC
h = [-0.25 -0.25 -0.25; -0.25 1 0; 0 0 0];
dc = reshape(dc, size(table) / 8);
dcdindex = floor(filter2(h, dc)) + 256;

% Get NZ and VV sequence of AC component.
[nz vv] = runlength(ac(:)', 63);

% Encode DC, ACNZ and ACVV components with the huffman code.
codeddcd = encodemessage(dcdindex(:)', dcdcode);
codednz = encodemessage(nz(:)' + 2, nzcode);
codedvv = encodemessage(vv(:)' + 129, vvcode);


% Header with codedacnz and codedacvv length
nbits = ceil(log2(63 * prod(size(table)) / 64));
header = strcat(dec2bin(length(nz), nbits), dec2bin(length(vv), nbits));

% After entropy coding, the final bitstream is constructed by concatenating
% the Huffman coding of DC and AC coefficients. In the example we have also
% included the header information with the input image size:

% Final bitstream
bitstream = strcat(header, codeddcd, codednz, codedvv);

% In the case of the _Table Tennis_ image, the resulting size in bits of the final
% bitstream (without header information) corresponds to:

final_bits = size(bitstream,2),

% which yields the following bits per pixel used in the representation of
% the grey image:

bitspixel = final_bits / (size(table,1)*size(table,2)),

% If we compare the obtained bits per pixel with the case of using 8 bits
% per pixel (when there is no compression), the final compression ratio
% achieved by our JPEG oriented encoder for the Table Tennis image is:

compression_ratio = 8 / bitspixel

% In order to compare the final compressed image with the original one, the
% produced bitstream must be decoded. The following code shows how the
% bitstream can be decoded back into quantized DCT coefficients.

% Extract header from codedmessage and calculate number of ACNZ and ACVV to decode.
zzsize = size(table) .* [1/8 8];
dcsize = [zzsize(1, 1) ceil(zzsize(1, 2) / 64)];
ndc = prod(dcsize); nac = ndc * 63;
nbits = ceil(log2(nac));
nzlen = bin2dec(bitstream(1 : nbits));
vvlen = bin2dec(bitstream(nbits + 1 : 2 * nbits));
bitstream2 = bitstream(2 * nbits + 1 : end);

load('decoders/im/dcddecoder.mat');
[dcdindex bitstream2] = decodemessage(bitstream2, dcddecoder, ndc);
load('decoders/im/nzdecoder.mat');
[nzindex bitstream2] = decodemessage(bitstream2, nzdecoder, nzlen);
load('decoders/im/vvdecoder.mat');
[vvindex bitstream2] = decodemessage(bitstream2, vvdecoder, vvlen);

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

% reconstruct the image
fun = @(x) idct2(iquantizedct8(izigzag8(x),Q));
tabler = uint8(blkproc(zzdctq, [1 64], fun));

% The PSNR is a good measure of the visual quality of the resulting
% compressed image (see Section 8.1.5). In the case of the Table Tennis
% image that we have used in all our examples, the final PSNR (for k = 1)
% is:

mypsnr(table,tabler),


%% Still image coding

% For completeness, the MATLAB function stillimagegrayencode implements the
% total encoding of a single image following the concepts explained in this
% Chapter. It accepts the matrix of the original image and the factor k to
% multiply the quantization table. The function returns the reconstructed
% image, the bitstream with the encoded DCT coefficients and the peak
% signal-to-noise ratio (PSNR) of the compressed image.

% Fig. 8.30 presents the result of applying this function to the Foreman
% CIF and Tennis Table QCIF original images previously used.

table = imread('seqs/table/table_000_g.bmp');
for k=[1 3 5],
    [Ir,bits,p] = stillimagegrayencode(table,k);
    figure;imshow(Ir);
    tb = size(bits,2) / numel(table);
    disp(['k = ' num2str(k) ' psnr = ' num2str(p) ' bits/pixel = ' num2str(tb) ' compression ratio = ' num2str(8/tb)]);
end;

foreman = imread('seqs/foreman/foreman_cif_000_g.bmp');
for k=[1 3 5],
    [Ir,bits,p] = stillimagegrayencode(foreman,k);
    figure;imshow(Ir);
    tb = size(bits,2) / numel(foreman);
    disp(['k = ' num2str(k) '; psnr = ' num2str(p) ' dB bits/pixel = ' num2str(tb) ' compression ratio = ' num2str(8/tb)]);
end;

% For color images, the same concepts apply to each of the three color
% channels (in the case of JPEG standard YCbCr color space is commonly used
% and the chrominance channels are decimated by factor of 2 in horizontal
% and vertical directions). The MATLAB function stillimageencode implements
% the encoder for color images. Fig. 8.31 shows the resulting compressed
% images for the same Table Tennis (QCIF) and Foreman images of previous
% examples (in QCIF and CIF formats), respectively:

table = imread('seqs/table/table_000.bmp');
for k=[1 3 5],
    [Ir,bits,p] = stillimageencode(table,k);
    figure;imshow(Ir);
    tb = size(bits,2) / size(table,1) / size(table,2);
    disp(['k = ' num2str(k) ' psnr = ' num2str(p) ' bits/pixel = ' num2str(tb) ' compression ratio = ' num2str(8*3/tb)]);
end;

table = imread('seqs/foreman/foreman_000.bmp');
for k=[1 3 5],
    [Ir,bits,p] = stillimageencode(table,k);
    figure;imshow(Ir);
    tb = size(bits,2) / size(table,1) / size(table,2);
    disp(['k = ' num2str(k) ' psnr = ' num2str(p) ' bits/pixel = ' num2str(tb) ' compression ratio = ' num2str(8*3/tb)]);
end;

table = imread('seqs/foreman/foreman_cif_000.bmp');
for k=[1 3 5],
    [Ir,bits,p] = stillimageencode(table,k);
    figure;imshow(Ir);
    tb = size(bits,2) / size(table,1) / size(table,2);
    disp(['k = ' num2str(k) ' psnr = ' num2str(p) ' bits/pixel = ' num2str(tb) ' compression ratio = ' num2str(8*3/tb)]);
end;
