%% Chapter 3 - How is sound processed in an MP3 player?
% This is a companion file to the book "Applied Signal Processing", 
% by T.Dutoit and F. Marques, Springer 2008.
% 
% It is supposed to be run cell-by-cell, using the cell mode of 
% MATLAB 6 and later versions. Search for "what are cells" in the 
% product help of MATLAB to know how to use them.
% 
% This file uses the SIGNAL_PROCESSING toolbox of MATLAB.

[input_signal,Fs]=audioread('violin.wav');
N = length(input_signal)
slen = N/Fs

[y,Fs]=audioread('01 David Bowie - Blackstar.wav');
cut = y(1:slen*Fs,:);
mono = sum(cut,2)/2;
 

%% 3. 32-Channel Pseudo-QMF filter bank
% We now build a 32-channel PQMF filter bank, as implemented in the MPEG-1
% Layer-I norm, and check its perfect reconstruction capability.
%
% *MATLAB function involved:*
% 
% * |hn = PQMF32_prototype| returns in |hn| the impulse response of the prototype
% low-pass symmetric filter of length 512 for building a 32-channel PQMF
% filter bank. This filter is used in the MPEG-1 Layer-I coder. Its
% normalized bandpass is 1/64 Hz and it satisfies the PR condition.  

% Load the prototype lowpass filter
hn=PQMF32_prototype;

%%
% Build 32 cosine modulated filters centered on normalized frequencies
% Fi=(2*i+1)/64 *1/2
PQMF32_Gfilters = zeros(32, 512);
for i = 0:31
    t2 = ((2*i+1)*pi/(2*32))*((0:511)+16);
    PQMF32_Gfilters(i+1,:) = hn.*cos(t2);
end
PQMF32_Hfilters=fliplr(PQMF32_Gfilters);

for i = 1:32
    [H,W]=freqz(PQMF32_Hfilters(i,:),1,512,44100);
    plot(W,20*log10(abs(H))); hold on
end
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
set(gca,'xlim',[0 22050]);
hold off;

%%
% Clearly, PQMF filters are not ideal band-pass filters with normalized
% bandwidth=1/64. But they almost satisfy the "pseudo" feature, which
% assumes that the filter response in a given band only overlaps with its 2
% neighboring bands. The total normalized bandwidth of each filter is 
% indeed about  1/32. Notice the filter gain is 15 dB, i.e. 20*log10(32)/2.
% As a result, passing twice through the filter produces a gain of 32,
% which compensates for the decimation by 32. In the next lines, it is thus
% not necessary to multiply upsampled sub-band signals by 32.
%
% Let us now check the output of the PQMF filter bank when fed with 2
% seconds of violin monophonic signal sampled at 44.100 Hz.

%[mono,Fs]=audioread('01 David Bowie - Blackstar.wav');
clf;
specgram(mono,1024,Fs,256);
soundsc(mono,Fs);

%%
output_signal=zeros(size(mono));

for i=1:32

    Hi_output=filter(PQMF32_Hfilters(i,:),1,mono);
    subband_i=Hi_output(1:32:end);
    upsampled_i(1:32:32*length(subband_i))=subband_i;
    % synthesis filters are the symmetric of the analysis filters, which
    % ensures linear phase in each sub-band
    Gi_output=filter(PQMF32_Gfilters(i,:),1,upsampled_i');
    output_signal=output_signal+Gi_output;
    
    % Screen output
    fprintf('processing sub-band %3d\n', i)

    if i==3 
        G3_output=Gi_output;
    end;

end;

%%
% As revealed by listening sub-band 3, isolated sub-band signals
% are very much aliased, because each  PQMF filter is not ideal. 

clf;
specgram(G3_output,1024,Fs,256);
soundsc(G3_output,Fs);

%%
% The PQMF filter bank makes sure aliasing in adjacent bands cancels itself
% when sub-bands are added. 

specgram(output_signal,1024,Fs,256);
soundsc(output_signal,Fs);

%%
% The power of the reconstruction error is about 85 dB below that of the
% signal. Notice the ouput is delayed by 511 samples (since H and G filters
% have a delay of 511/2 samples).
%
% *MATLAB function involved:*
% 
% * |snr = snr(signal,signal_plus_noise, max_shift,showplot)| returns the
% signal-to-noise ratio computed from the input signals. |Max_shift| gives
% the maximum time-shift (in samples) between |signal| and
% |signal_plus_noise|. The actual time-shift (obtained from the maximum of
% the cross-correlation between both signals) is taken into account to
% estimate the noise. If |showplot| is specified, then the signal,
% signal_plus_noise, and error are plotted, and the SNR is printed on the
% plot.

error=output_signal(512:end)-mono(1:end-511);

[signal_psd,w]=periodogram(mono(11001:12024),hamming(1024));
[error_psd,w]=periodogram(error(11001:12024),hamming(1024));
plot(w/pi*22050,10*log10(signal_psd)); hold on;
plot(w/pi*22050,10*log10(error_psd),'r','linewidth',2); hold off;
set(gca,'xlim',[0 22050]);
legend('Signal PSD', 'Error PSD');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');

snr_PQMF=snr(mono(1:end-511),output_signal(512:end),0)

%% 4. Filter banks and Lapped Transforms
% If the length of the filters used in M-channel sub-band filters were
% exactly M samples, one could easily see the operations performed
% simultaneously by the M analysis filters as a linear transform of
% successive non-overlapping M sample frames of the input signal.
%
% A 4-sample DFT, for instance, can implement a 4-channel filter bank
% whose subband filters are the time-reversed of the lines of the 4x4 DFT
% matrix. Applying it to a chirp is straightforward.

Fs=8000;
mono=chirp((1:4*Fs)/Fs,0,4,4000);

for i=1:length(mono)/4

    % creating a column vector with 4 samples
    input_frame=mono(4*(i-1)+1:4*(i-1)+4)';
    
    % producing one sample in each downsampled sub-band, i.e. band-pass
    % filtering and downsampling all sub-bands in one operation.
    subbands_sample=fft(input_frame);
    
    % producing four samples of the filter bank output, i.e. upsampling, 
    % band-pass filtering all sub-bands, and summing them in one operation. 
    output_frame=ifft(subbands_sample);
    
    % storing the output column vector in the output signal
    output_signal(4*(i-1)+1:4*(i-1)+4)=output_frame';

end;

%soundsc(output_signal,Fs);

%%
% Since the underlying filters have complex coefficients, however, each
% sub-band signal is complex. What is more, this type of filter bank is not
% very frequency selective, as shown below. The frequency overlap between
% adjacent bands is about half the main lobe band-pass (as in the previous
% section on PQMF), but the side lobes are very high. This does not make it
% a good candidate for sub-band coding. 

tmp=[1 ; exp(-j*pi/2) ; exp(-j*pi) ; exp(-j*3*pi/2)];
DFT_matrix_4x4=vander(tmp); 

for i=1:4
    [H,W]=freqz(fliplr(DFT_matrix_4x4(i,:)),1,'whole');
    plot(W/pi,max(20*log10(abs(H)),-50)); hold on
end
xlabel('Normalized frequency (*pi rad/sample)'); ylabel('Magnitude (dB)');
hold off;

%%
% In general, the length of the impulse responses of the analysis and
% synthesis filters used in sub-band coders is higher than the number M of
% channels.  The filtering operations, however, can still be
% implemented as the multiplication of L-sample frames with LxM or
% MxL matrices. 
%
% For example, the 32-channel PQMF filter-bank introduced in the previous
% Section, in which the length L of the impulse response of each filter is
% 512 samples, can be efficiently implemented as follows (which is very much
% similar to the implementation of our previous DFT-based filter bank,
% with the addition of overlap).

% Build the PQMF H and G filters
hn=PQMF32_prototype;
PQMF32_Gfilters = zeros(32, 512);
for i = 0:31
    t2 = ((2*i+1)*pi/(2*32))*((0:511)+16);
    PQMF32_Gfilters(i+1,:) = hn.*cos(t2);
end

[mono,Fs]=audioread('violin.wav');

% Block-based sub-band filtering
input_frame=zeros(512,1);
output_signal=zeros(size(mono));

for i=1:(length(mono)-512+32)/32
    
    % Overlap input_frames (column vectors)
    input_frame=mono((i-1)*32+1:(i-1)*32+512);

    % Analysis filters and downsampling
    % Since PQMF H filters are the time-reversed G filters, we use the G
    % filters matrix to simulate analysis filtering
    subbands_frame_i = PQMF32_Gfilters*input_frame;
    
    % Synthesis filters
    output_frame = PQMF32_Gfilters'*subbands_frame_i;

    % Overlap output_frames (with delay of 511 samples)
    output_signal((i-1)*32+1:(i-1)*32+512)= ...
        output_signal((i-1)*32+1:(i-1)*32+512)+output_frame;
    
end

%%
% Obviously we get the same results as before, and the overall SNR is
% unchanged.

error=output_signal-mono;

[signal_psd,w]=periodogram(mono(11001:12024),hamming(1024));
[error_psd,w]=periodogram(error(11001:12024),hamming(1024));
plot(w/pi*22050,10*log10(signal_psd)); hold on;
plot(w/pi*22050,10*log10(error_psd),'r','linewidth',2); hold off;
set(gca,'xlim',[0 22050]);
legend('Signal PSD', 'Error PSD');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');

snr_lapped=snr(mono(512:end-512),...
                     output_signal(512:end-512),0)

%% 5. Perceptual audio coding
% The sub-band filtering process developed in the previous Sections
% transforms the original stream of samples at sampling frequency Fs into
% 32 parallel sub-bands sampled at Fs/32. It does not by itself produce
% compression.
%
% Quantizing sub-band samples uniformly does not allow much transparency at
% low bit rates, as shown in the 4-bits per sub-band sample trial below
% (compression factor = 4). We use a mid-thread quantizer here, so as to
% encode low level signals to 0.

% Build the PQMF H and G filters
hn=PQMF32_prototype;
PQMF32_Gfilters = zeros(32, 512);
for i = 0:31
    t2 = ((2*i+1)*pi/(2*32))*((0:511)+16);
    PQMF32_Gfilters(i+1,:) = hn.*cos(t2);
end

[mono,Fs]=audioread('violin.wav');

% Block-based sub-band analysis filtering
input_frame=zeros(512,1);
output_signal=zeros(size(mono));

n_frames=(length(mono)-512+32)/32;
for i=1:n_frames
     
     % Overlap input_frames (column vectors)
     input_frame=mono((i-1)*32+1:(i-1)*32+512);
 
     % Analysis filters and downsampling
     % NB: we put sub-band signals in columns
     subbands(i,:) = (PQMF32_Gfilters*input_frame)';
 
     % Uniform quantization on 4 bits, using a mid-thread quantizer in
     % [-1,+1] 
     n_bits = 4;
     alpha = 2^(n_bits-1);
     quantized_subbands(i,:) = (floor(alpha*subbands(i,:)+0.5))/alpha; % mid-thread

     % the |uencode| and |udecode| functions provided by MATLAB do not
     % properly implement a mid-thread quantizer. Using them here (by
     % uncommenting the next two lines results in tonal quantization noise.
     % See Appendix 1 at the end of this script. 
     % codes = uencode(subbands(i,:),4,1);
     % quantized_subbands(i,:) = udecode(codes,4);
 
     % Synthesis filters
     output_frame = PQMF32_Gfilters'*quantized_subbands(i,:)';
 
     % Overlap output_frames (with delay of 511 samples)
     output_signal((i-1)*32+1:(i-1)*32+512)= ...
         output_signal((i-1)*32+1:(i-1)*32+512)+output_frame;
     
end

specgram(output_signal,1024,Fs,256);
soundsc(output_signal,Fs);

%%
% The resulting output signal is degraded. It exhibits a strong high
% frequency tonal noise. This is typical of sub-band coding, in which
% quantization errors are produced at Fs/32, i.e. 1378 Hz.
% The SNR falls down to 10.3 dB.

% NB: no delay compensation required, as the first sub-band samples
% produced by the lapped transform correspond to the 512th original
% samples.

error=output_signal-mono;

[signal_psd,w]=periodogram(mono(11001:12024),...
    hamming(1024));
[error_psd,w]=periodogram(error(11001:12024),...
    hamming(1024));
plot(w/pi*22050,10*log10(signal_psd)); hold on;
plot(w/pi*22050,10*log10(error_psd),'r','linewidth',2); hold off;
set(gca,'xlim',[0 22050]);
legend('Signal PSD', 'Error PSD');
xlabel('Frequency (Hz)'); 
ylabel('Magnitude (dB)');

snr_4bits=snr(mono(512:end-512),output_signal(512:end-512),0)

%%
% One can see that the fixed [-1,+1] quantizer range does not adequately
% account for the  variation of sub-band signal level across sub-bands, 
% as well as in time.  

plot(quantized_subbands(100:200,2)); hold on; 
plot(subbands(100:200,2),'--r'); hold off;
legend('Sub-band #2, quantized', 'Sub-band#2, original');
xlabel('Time (samples at Fs/32)'); ylabel('Amplitude');

%%
% An obvious means of enhancing its quality is therefore to apply a
% scale factor to each sub-band quantizer. As in the MPEG-1 Layer-I coder,
% we compute a new scale factor every 12 sub-band sample (i.e. every 
% 32x12 = 384 sample at the original sample rate).
% Quantization errors are much reduced (here in sub-band #2).

% Adaptive quantization per blocks of 12 frames
n_frames=fix(n_frames/12)*12;
for k=1:12:n_frames
    
    % Computing scale factors in each 12 samples sub-band chunk
    [scale_factors,tmp]=max(abs(subbands(k:k+11,:)));
    
    % Adaptive uniform quantization on 4 bits, using a mid-thread quantizer in
    % [-Max,+Max] 
    for j=1:32 % for each sub-band
 
        n_bits = 4;
        alpha = 2^(n_bits-1)/scale_factors(j);
        quantized_subbands(k:k+11,j) = ...
            (floor(alpha*subbands(k:k+11,j)+0.5))/alpha; % mid-thread

    end;

end;

plot(quantized_subbands(100:200,2)); hold on; 
plot(subbands(100:200,2),'--r'); hold off;
legend('Sub-band #2, quantized', 'Sub-band#2, original');
xlabel('Time (samples at Fs/32)'); ylabel('Amplitude');

%%
% The resulting signal is of much higher quality.

% Signal synthesis
output_signal=zeros(size(mono));
for i=1:n_frames

    % Synthesis filters
    output_frame = PQMF32_Gfilters'*quantized_subbands(i,:)';

    % Overlap output_frames (with delay of 511 samples)
    output_signal((i-1)*32+1:(i-1)*32+512)= ...
        output_signal((i-1)*32+1:(i-1)*32+512)+output_frame;
    
end

specgram(output_signal,1024,Fs,256);
soundsc(output_signal,Fs);

%%
% The overall SNR has increased to 25 dB.

error=output_signal-mono;

[signal_psd,w]=periodogram(mono(11001:12024),...
     hamming(1024));
[error_psd,w]=periodogram(error(11001:12024),...
     hamming(1024));
plot(w/pi*22050,10*log10(signal_psd)); hold on;
plot(w/pi*22050,10*log10(error_psd),'r','linewidth',2); hold off;
set(gca,'xlim',[0 22050]);
legend('Signal PSD', 'Error PSD');
xlabel('Frequency (Hz)'); 
ylabel('Magnitude (dB)');

snr_4bits_scaled=snr(mono(512:end-512),...
    output_signal(512:end-512),0)

%%
% The ultimate refinement, which is by far the most effective and results
% from years of audio research, consists in accepting more quantization
% noise (by allocating less bits, and thereby accepting a higher SNR) in
% frequency bands where it will not be heard, and using these extra bits
% for more perceptually prominent bands. The required perceptual
% information is provided by a psycho-acoustical model.
%
% For any 512-sample frame taken from the input signal, the MPEG-1 Layer-I
% psycho-acoustical model computes a global masking threshold, obtained by
% first detecting prominent tonal and noise maskers, separately, and
% combining their individual thresholds. The maximum of this global
% threshold and the absolute auditory threshold is then taken as the final
% threshold. 
%
% *MATLAB function involved:*
% 
% * |function [SMR, min_threshold_subband, masking_threshold] = ...
%    MPEG1_psycho_acoustic_model1(frame)|
% Computes the masking threshold (in dB) corresponding to psycho-acoustic
% model #1 used in MPEG-1 Audio (cf  ISO/CEI  norm 11172-3:1993 (F), pp.
% 122-128). 
% Input |frame| length should be 512 samples. 
% |min_threshold_subband| returns the minimun of |masking threshold| in each
% of the 32 sub-bands. 
% |SMR| returns 27 signal-to-mask ratios (in dB);
% |SMR|(28-32) are not used. 

frame=mono(11001:11512);
[SMR, min_threshold,frame_psd_dBSPL]= ...
    MPEG1_psycho_acoustic_model1(frame);
% NB: the power levels returned by this function assume that a full-scale 
% signal (in [-1,+1]) corresponds to 96 DB SPL

f = (0:255)/512*44100;
auditory_threshold_dB = 3.64*((f/1000).^-0.8) - ...
    6.5*exp(-0.6.*((f/1000)-3.3).^2) + 0.001*((f/1000).^4);
plot(f, frame_psd_dBSPL, f, min_threshold,'.r', ...
    f, auditory_threshold_dB, '-.k');
hold off;
axis([0 22050 -20 100]);
legend('Signal PSD', 'Min. threshold per sub-band','Absolute threshold');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');

%%
% Signal-to-Mask Ratios (SMR) are computed for each band, in a very
% conservative way (as the ratio of the maximum of the signal PSD to the
% minimum of the masking threshold in each band).  
% Bit allocation is performed by an iterative algorithm which gives
% priority to sub-bands with higher SMR. 
% The resulting SNR in each sub-band should be greater or equal to
% the SMR, so as to push the noise level below the masking threshold.
%
% *MATLAB function involved:*
% 
% * |function [N_bits,SNR] = MPEG1_bit_allocation(SMR, bit_rate)|
% Implements a simplified bit allocation greedy algorithm.
% |SMR| is the signal-to-mask ratios in each sub-band,
% as defined by the MPEG1 psycho-acoustic model. 
% |bit_rate| is in kbits/s.
% |N_bits| is the number is bits in each sub_band.
% |SNR| is the maximum SNR in each sub-band after quantization, i.e. the
% SNR assuming each sub-band contains a full-range sinusoid.
% NB: N_bits and SNR are set to zero for sub-bands 28 to 32.

% Allocating bits for a target bit rate of 192 kbits/s (compression
% ratio = 4)
[N_bits, SNR] = MPEG1_bit_allocation(SMR, 192000);

stairs((0:32)/32*22050,[SMR SMR(32)]);
hold on;
stairs((0:32)/32*22050,[SNR SNR(32)],'--');
axis([0 22050 -20 100]);
legend('SMR', 'SNR');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');

%%
% Let us test this on the complete signal, adding perceptual bit allocation
% to adaptive uniform quantization. 
% Notice we simplify the quantization
% scheme here, compared to MPEG-1, by considering quantization with any
% number of bits in [0,16].  

% Adaptive quantization per blocks of 12 frames
n_frames=fix(n_frames/12)*12;
for k=1:12:n_frames
    
    % Computing scale factors in each 12 samples sub-band chunk
    [scale_factors,tmp]=max(abs(subbands(k:k+11,:)));
    
    % Computing SMRs
    frame=mono(176+(k-1)*32:176+(k-1)*32+511);
    % NB: the input frame for the psycho-acoustic model is delayed by 176
    % samples, as it should be centered in the middle of the local block of
    % 12 sub-band samples, and the first sub-band sample corresponds to
    % original sample 256 (actually, to sample 512, but the analysis filter
    % introduces a delay of 256 samples). So the center of the first frame should
    % be on original sample 256+11*32/2=432, and the beginning of this frame
    % falls on sample 11*32/2=176.
    SMR = MPEG1_psycho_acoustic_model1(frame);

   % Allocating bits for a target bit rate of 192 kbits/s (compression
   % ratio = 4)
    N_bits = MPEG1_bit_allocation(SMR, 192000);

    % Adaptive perceptual uniform quantization, using a mid-thread 
    % quantizer in [-Max,+Max] 
    for j=1:32 % for each sub-band

        if N_bits(j)~=0
            alpha = 2^(N_bits(j)-1)/scale_factors(j);
            quantized_subbands(k:k+11,j) = ...
                (floor(alpha*subbands(k:k+11,j)+0.5))/alpha; % mid-thread

%             codes= uencode(subbands(k:k+11,j),N_bits(j),...
%                 scale_factors(j),'signed');
%             quantized_subbands(k:k+11,j) = ... 
%                 udecode(codes,N_bits(j),scale_factors(j));
        else
            quantized_subbands(k:k+11,j) = 0;
        end;

    end;

    % Screen output
    fprintf('processing frame %3d/%3d\n', k,n_frames)

end;

% Signal synthesis
output_signal=zeros(size(mono));
for i=1:n_frames
    % Synthesis filters
    output_frame = PQMF32_Gfilters'*quantized_subbands(i,:)';
    % Overlap output_frames (with delay of 511 samples)
    output_signal((i-1)*32+1:(i-1)*32+512)= ...
        output_signal((i-1)*32+1:(i-1)*32+512)+output_frame;
end

clf;
specgram(output_signal,1024,Fs,256);
soundsc(output_signal,Fs);

%%
% Quantization noise is now very small in some prominent sub-bands, like
% sub-band #2. 

clf;
plot(subbands(100:200,2),'--r'); hold on;
plot(quantized_subbands(100:200,2)); hold off; 
legend('Sub-band #2, original','Sub-band#20, quantized');
xlabel('Time (samples at Fs/32)'); ylabel('Amplitude');

%%
% It is also more important in other sub-bands, like sub-band #20. 

plot(subbands(100:200,20),'--r'); hold on;
plot(quantized_subbands(100:200,20)); hold off; 
legend('Sub-band #20, original','Sub-band#20, quantized');
xlabel('Time (samples at Fs/32)'); ylabel('Amplitude');

%%
% The overall SNR has increased to 36.2 dB.
% The perceptual SNR is actually much higher, since most of the noise
% cannot be heard. 

error=output_signal(1:end)-mono(1:end);

[signal_psd,w]=periodogram(mono(11001:12024),...
    hamming(1024));
[error_psd,w]=periodogram(error(11001:12024),...
    hamming(1024));
plot(w/pi*22050,10*log10(signal_psd)); hold on;
plot(w/pi*22050,10*log10(error_psd),'r','linewidth',2); hold off;
legend('Signal PSD', 'Error PSD');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');

snr_scaled_perceptual=snr(mono(512:end-512),...
    output_signal(512:end-512),0)

%%
% The price to pay is that the scale factors and the number of bits per
% sub-band must now be stored, every 12 frames (i.e. every 12*32 sub-band
% samples). 
% In the MPEG-1 norm, scale factors expressed in dB and quantized
% on 6 bits each. This comes from the fact that the ear perceives loudness
% as the log of the energy, and has about 96 dB of hearing dynamics with a
% sensitivity threshold of about 1 dB. With 6 bits, the quantization step
% is 96/64 dB and the error lies in [-96/62/2, +96/64/2], which is below
% the 1 dB threshold.
% Assuming 4 bits are required for this information in each of the
% 27 first sub-bands (sub-bands 28 to 32 are ignored my MPEG1), this leads
% to 27*10=270 bits every 12 frames. 
% (In practice, MPEG-1 does not allow all integer values in [1,16] for the
% number of bits in each band, which makes it possible to quantize it to
% less then 4 bits. The total number of bits used for bit allocation is
% then reduced to 88.) 
% It is easy to obtain the number of bits used for the last block of 12
% frames in our audio test, and the related bit-rate.

bits_per_block=sum(N_bits)*12+270
bit_rate=bits_per_block*44100/384

%% Appendix 1: Uniform quantizers
% In Section 5 we have used a mid-thread quantizer. It was tempting to use
% the |uencode| and |udecode| fucntions provided by MATLAB, but a quick
% examination of the quantization charactertic of these functions shows
% they do not implement a real mid-thread quantizer.

x = -1:0.01:1;
n_bits=3;
codes = uencode(x,n_bits,1);
y1 = udecode(codes,n_bits,1);
alpha = 2^(n_bits-1);
y2 = (floor(alpha*x)+0.5)/alpha; % mid-rise
y3 = (floor(alpha*x+0.5))/alpha; % mid-thread

plot(x,y1); hold on
plot(x,y2,'r')
plot(x,y3,'g')
plot([0 0], [-1 1])
plot([-1 1], [0 0])
axis([-1 1 -1 1])
legend('uen(de)code', 'mid-rise','mid-thread');