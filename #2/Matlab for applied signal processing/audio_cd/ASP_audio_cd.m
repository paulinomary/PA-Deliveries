%% Chapter 2 - How are bits played back from an audio CD?
% This is a companion file to the book "Applied Signal Processing", 
% by T.Dutoit and F. Marques, Springer 2008.
% 
% It is supposed to be run cell-by-cell, using the cell mode of 
% MATLAB 6 and later versions. Search for "what are cells" in the 
% product help of MATLAB to know how to use them.
% 
% This file uses the SIGNAL_PROCESSING toolbox of MATLAB.
 
%%
% In this chapter we will see that, contrary to what might be expected,
% many audio digital-to-analog converters (DACs), including those used in
% CD and MP3 players, first requantize the 16-bit stream into a one-bit stream
% with very high samping frequency, using a signal processing concept known
% as the delta-sigma modulation, and then convert the resulting binary
% signal back to an audio waveform. 
% 
% We will first revise the basics of uniform quantization (Section 1), including
% the important part played by dithering (Section 2). We will then compare
% the internals of a conventional DAC (Section 3) to those of more advanced
% DACs, using oversampling (Section 4), noise-shaping (Section 5), and
% delta-sigma modulation (Section 6).
%
% Copyright T. Dutoit, R. Schreier (2007)

set(0,'defaultFigureColor','w');

%% 1. Uniform quantization
% In order to perform D/A conversion, we first need to create a PCM signal.
% We will first check the main results of uniform quantization theory,
% which we will use later.
%
% Let us first generate 8000 samples of uniform white noise with zero mean,
% which we will assume to have been sampled at Fs=8000 Hz. For convenience,
% we set its variance to 1 (0 dB) by imposing its range to [-peak,+peak],
% with peak=sqrt(12)/2 ~ 1.73. This is confirmed on the power spectral
% density of the signal. 

signal_std=1;
peak=sqrt(12*signal_std)/2;
signal=(rand(1,8000)-0.5)*2*peak;

pwelch(signal,[],[],[],2);
soundsc(signal,8000); 

%%
% (For comments on  MATLAB's |pwelch| function and the way we use it here,
% refer to Appendix I at the end of this document.) 

%% 
% Quantizing this signal uniformly on N bits is easy. Let us choose N=3 for
% convenience. We use mid-rise quantization (hence the |+q/2|), best suited
% for an even number of quantization steps.

N=3;
quantizer_saturation=peak;
q=(2*quantizer_saturation)/2^N;
signal_quantized=floor(signal/q)*q+q/2; % mid-rise

plot(signal(1:50)); hold on; 
stem(signal_quantized(1:50),'.');hold off
xlabel('Time (samples)'); ylabel('Amplitude');
legend('signal','quantized signal');

%%
% As expected, the quantization error (or noise) is in the range
% [-q/2,q/2].

error=signal-signal_quantized;

plot(signal(1:50),'-.'); hold on; 
plot(error(1:50));hold off
xlabel('Time (samples)'); ylabel('Amplitude');
legend('signal','quantization error');

%%
hist(error);

%%
% The quantization noise is also white, as revealed by power spectral
% density estimation. The variance of the noise can thus
% be read (in dB) on the psd plot. From the load factor (i.e., the
% saturation value of the quantizer, normalized by the standard deviation 
% of the signal), we can compute the theoretical SNR due to quantization. 
% It matches to the SNR computed from the samples, and corresponds to 
% the variance of the noise (with opposite sign) since the variance 
% of the signal is 0dB. 
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

pwelch(error,[],[],[],2);
load_factor=quantizer_saturation/signal_std;
snr_theory=6.02*N + 4.77 - 20*log10(load_factor)
snr_measured=snr(signal,signal_quantized,0)

%% 2. Dithering
% When the amplitude of the input signal approaches the quantization step,
% the quantization error becomes correlated with the signal. If the signal
% itself is not random, the quantization error can then be heard as
% non-linear audio distortion (rather than as additive noise).
%
% Let us check this by first creating a sine wave with fundamental
% frequency  f0=200 Hz, sampled at Fs=8000 Hz, and set its variance to 1 (0
% dB) by imposing its peak to sqrt(2). 

signal_std=1;
peak=sqrt(signal_std*2);
signal=peak*sin(2*pi*200*(0:1/8000:2));
var_signal_dB=10*log10(var(signal))
soundsc(signal,8000);

%%
% When we quantize it uniformly to 3 bits in the range [-2,+2], we notice
% that the quantization error does not look like real noise.
 
N=3;
quantizer_saturation=2;
q=2*quantizer_saturation/2^N;
signal_quantized=floor(signal/q)*q+q/2; % quantization
error=signal-signal_quantized;

plot(signal(1:50),'-.'); hold on; 
plot(error(1:50));hold off
xlabel('Time (samples)'); ylabel('Amplitude');
legend('signal','quantization error');

%% 
% Harmonics of the original frequency appear in the spectrum of the
% quantized signal and sound very unpleasant. The SNR is not directly
% readable on the psd plot, because of these harmonics 
%
% (NB: real frequencies = annouced frequencies * Fs/2 => the original 200
% Hz signal is shown at F=0.05) 

snr_quantized=snr(signal,signal_quantized,0)
pwelch(signal_quantized,[],[],[],2);
soundsc(signal_quantized,8000);

%%
% In such a case it may be interesting to whiten the quantization error by
% adding real noise to it. This operation is termed as dithering.
%
% Let us add a dithering noise with triangular pdf in the range [-q,+q],
% i.e. two times the quantization step. Such a noise is easily obtained by
% adding two uniform white noises in the range [-q/2,q/2].

dither=(rand(size(signal))-0.5)*q+(rand(size(signal))-0.5)*q;
hist(dither);
signal_dithered=signal+dither;

%%
% The resulting quantization error looks indeed more noise-like.
% (Notice we do not account for quantizer saturation, which 
% may occur when dither is added to a signal which approaches 
% the saturation levels of the quantizer. This is note the case 
% in our example.)

signal_dithered_quantized=floor(signal_dithered/q)*q+q/2;
error_dithered_quantized=signal-signal_dithered_quantized;

plot(signal(1:50),'-.'); hold on; 
plot(error_dithered_quantized(1:50));hold off
xlabel('Time (samples)'); ylabel('Amplitude');
legend('signal','quantization error, with dithering');

%% 
% The power spectral density of the quantized signal now appears as white
% noise added to the initial sine wave.

pwelch(signal_dithered_quantized,[],[],[],2);

%%
% Dithering clearly degrades the SNR (by about 4.8dB), but results in
% perceptually more acceptable quantization error.

load_factor=quantizer_saturation/signal_std;
snr_theory=6.02*N + 4.77 - 20*log10(load_factor)
snr_dithered_quantized=snr(signal,signal_dithered_quantized,0)

soundsc(signal_dithered_quantized,8000);

%% 3. Conventional DAC
% A conventional DAC first creates a staircase signal from the sequence of
% samples, by converting each sample into an analog voltage and holding
% this voltage for 1/Fs seconds. This is called zero-order (analog)
% interpolation. This operation is critical: the higher the number of bits
% in the PCM code to convert, the higher the precision required for
% creating the staircase signal!
%
% We will demonstrate this on a sine wave, with N=6 bits. 
%
% *MATLAB function involved:*
% 
% * |signal_quantized = uquantize(signal,N,saturation)| quantizes |signal|
% uniformly to N bits in the range [-saturation, +saturation], using a mid-rise
% quantizer and triangular dither with a range of twice the quantization
% step.

signal_std=1;
peak=sqrt(signal_std*2);
signal=peak*sin(2*pi*200*(0:1/8000:2));

N=3;
quantizer_saturation=2;
signal_quantized=uquantize(signal,N,quantizer_saturation);
snr_quantized=snr(signal,signal_quantized,0)

%%
% We can SIMULATE the analog zero-order interpolation in the digital domain
% by working with a much higher sampling frequency than Fs (say, Fs'=10 Fs,
% i.e. 80 kHz). Zero-order interpolation is then equivalent to inserting 9
% zeros in between each sample, and convolving the resulting signal
% with a sequence of 10 unity samples. (This is actually a very rudimentary
% form of numerical interpolation).

signal_pulses=zeros(1,10*length(signal_quantized));
signal_pulses(1:10:10*length(signal))=signal_quantized;

pwelch(signal_pulses,[],[],[],20);

%%
hold_impresp=ones(1,10);
signal_staircase=conv(signal_pulses,hold_impresp);

plot(1:10:391, signal_quantized(1:40),'o'); hold on;
plot(signal_staircase(1:400)); hold off;
xlabel('Time (samples)'); ylabel('Amplitude');
legend('samples','after zero hold');

%%
% Notice that zero-order interpolation acts as a lowpass filter, performing
% a first attenuation of the spectral images of the signal at integer
% multiples of Fs.

freqz(hold_impresp,1, 256, 80000);

%%
% This directly affects the spectrum of the resulting "analog" signal. One
% can clearly see spectral images, weighted by the effect of the
% interpolation.

pwelch(signal_staircase,[],[],[],20);

%%
% The DA converter then feeds the resulting staircase signal to an analog
% low-pass smoothing filter, which removes the spectral images due to
% sampling. The passband of this filter is limited by the maximum frequency
% of the signal, Fm, and its stopband must not be greater than Fs-Fm. Let
% us assume that the PCM signal we have is a telephone signal, with Fm=3400
% Hz and Fs=8000 Hz, and perform the approximation of the filter with 0.1dB
% of ripple in the passband and -60 dB in the stopband. We use Chebyshev
% approximation, so as to keep the order of the filter low.

[order,Wn]=cheb1ord(2*pi*3400,2*pi*4600,0.1,60,'s');
[Num_LP,Den_LP]=cheby1(order,0.1,Wn,'s');
zplane(Num_LP,Den_LP);

%%
% The frequency response of this 12th order filter meets our requirements. 
% (Notice its phase is non-linear, but for audio applications this is not a
% problem). 

freqs(Num_LP,Den_LP);

%%
% We can now simulate analog filtering by convolving our highly oversampled
% staircase signal with a sampled version of the impulse response of the LP
% filter. 
% 
% *MATLAB function involved:*
% 
% |y = filters(N,D,x,Fs)| simulates analog filtering of the data in vector
% x with the filter N(s)/D(s) described by vectors N and D to create the
% filtered data y. Fs is the sampling frequency of the input (and output).
% Filtering is performed by convolving the input with an estimate of the
% impulse response of the filter, obtained by partial fraction expansion. 

analog_output=filters(Num_LP,Den_LP,signal_staircase,80000);

clf; plot(analog_output(1000:2000));
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% The noise level has been reduced by the analog LP filter,
% as revealed by spectral analysis. The final SNR is a bit higher than 
% what we had before D/A conversion, but this is due to the fact that the
% zero-order interpolation attenuates the noise at the upper edge of its
% passband.
%
% Notice, though, that the samples sent to the MATLAB's |soundsc| function
% are made available to your ears by ... another (real) DAC. One should
% therefore consider the final quality of this audio sample with care. 

pwelch(analog_output,[],[],[],20);
signal_sampled_at_10Fs=peak*sin(2*pi*200*(0:1/80000:2));
snr_analog =snr(signal_sampled_at_10Fs,analog_output,400) % 400=possible delay
soundsc(analog_output,80000);

%% 4. Oversampling DAC
% By oversampling the digital signal prior to sending it to the analog part
% of the DAC (which is responsible for creating a staircase analog signal
% before final analog low-pass filtering.), we can broaden the transition
% band of the analog smoothing filter, thereby strongly decreasing its
% order (and hardware complexity).
% 
% What is more, if we requantize the signal on N-1 bits (1 bit less=>6 more
% dBs of total noise power) after multiplying Fs by 4, only one quarter of
% the resulting power spectral density will contribute to quantization
% noise in the [0,Fs/2] band; the audible part of the requantization noise
% power will thus be 10log10(1/4)dB (i.e. 6 dB) lower than its total power.
% Hence, this N-1 requantization step will be felt as a new N-bit
% quantization. 
%
% In other words, dropping one bit from 4 times oversampled digital signals
% (or k bits for after oversampling by a ratio of k*4) does not do much
% harm to the signal, while it decreases the required precision of the
% hardware DAC  
%
% Let us check this on our sinusoidal test signal, quantized to 6 bits.

signal_std=1;
peak=sqrt(signal_std*2);
signal=peak*sin(2*pi*200*(0:1/8000:2));

N=6;
quantizer_saturation=2;
signal_quantized=uquantize(signal,N,quantizer_saturation);

pwelch(signal_quantized,[],[],[],2);
snr_quantized=snr(signal,signal_quantized,0)
soundsc(signal_quantized,8000);

%%
% We now oversample this quantized signal by a factor of 4, by first adding
% 3 zeros in between its samples. The resulting signal has a sampling
% frequency Fs'=32000 Hz. The amplification samples by 4 is required for
% the final sine wave to have the same peak-to-peak level as the original.

signal_pulses=zeros(1,4*length(signal_quantized));
signal_pulses(1:4:4*length(signal_quantized))=4*signal_quantized;

pwelch(signal_pulses,[],[],[],8);

%%
% Interpolation is performed by filtering the impulses with a quarter-band
% digital LP filter. Notice we use here a very sharp filter, which 
% requires a length of 300 coefficients. In oversampling DACs, LP filtering
% is performed much more crudely, for meeting low computational load
% constraints.

lp_fir=firpm(300,[0 .24 .26 1],[1 1 0 0]); % Linear phase quarter-band FIR filter
signal_interpolated=filter(lp_fir,1,signal_pulses);

plot(signal_interpolated(1000:1400));
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% The variances of the signal and of the quantization noise (hence, the
% SNR) have not changed.

pwelch(signal_interpolated,[],[],[],8);
signal_sampled_at_4Fs=peak*sin(2*pi*200*(0:1/32000:2));
snr_interpolated=snr(signal_sampled_at_4Fs,...
                            signal_interpolated,300)

%%
% Here comes the first positive effect of oversampling on SNR:
% (re)quantizing the interpolated signal with N-1 bit (i.e., omitting the
% least significant bit of the underlying PCM codes) does not affect the
% SNR. 
%
% Let us first perform 4-bit requantization of the interpolated signal.

signal_requantized=uquantize(signal_interpolated,4,...
                                        quantizer_saturation);

plot(signal_sampled_at_4Fs(2000:2400)); hold on;
plot(signal_requantized(2151:2550)); hold off;
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% The psd of the requantized signal shows noise with a variance of about
% 12dB below its initial value, but only 1/4th of that noise is in the band [0
% Fs/2]. Therefore the measured SNR does not reflect the actual SNR. We
% will see below that the actual SNR has only decreased by about 6 dB. 
% Physically speaking, this is perfectly understandable: the quantization
% noise samples produced at 4 Fs are averaged by lowpass filtering.

pwelch(signal_requantized,[],[],[],8);
snr_requantized=snr(signal_sampled_at_4Fs,...
                            signal_requantized,300)

%%
% Now we can again SIMULATE the analog part of D/A conversion (as in
% Section 3), using an "analog sampling frequency" of 3*32=96 kHz. 
%
% Oversampling has a second positive effect on this part: it allows for a
% much wider transition band for the lowpass smoothing filter : [3400Hz
% 28600Hz]. This results in a simpler filter, of order 4. 
%
% It even has a third interesting effect: the same lowpass filter can be
% used for a large range or values for Fs (as required for the sound card
% of a computer, for instance).

[order,Wn]=cheb1ord(2*pi*3400,2*pi*28600,0.1,60,'s'); 
[Num_LP,Den_LP]=cheby1(order,0.1,Wn,'s'); % relaxed filter
zplane(Num_LP,Den_LP);

%%
% Filtering is simulated here as in Section 3.

signal_pulses=zeros(1,3*length(signal_requantized));
signal_pulses(1:3:3*length(signal_requantized))=signal_requantized;
hold_impresp=ones(1,3);
signal_staircase=conv(signal_pulses,hold_impresp);
analog_output=filters(Num_LP,Den_LP,signal_staircase,96000);

plot(analog_output(2000:3300));
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% As a result of the analog filtering, the apparent SNR is now about 7 dB
% lower than for our initial 6-bit quantization, i.e. only 1dB less than
% 5-bit quantization at 8 kHz. The lost dB comes from dithering, and from
% the fact that the LP filter specifications are loose.

pwelch(analog_output,[],[],[],24);
signal_sampled_at_12Fs=peak*sin(2*pi*200*(0:1/96000:2));
snr_analog=snr(signal_sampled_at_12Fs,analog_output,4000)
soundsc(analog_output,96000);

%%
% Summarizing: a 6-bit DAC operating at Fs leads to the same SNR as 
% a 5-bit DAC operating at the 4x oversampling rate 4*Fs. 

%% 5. Oversampling and noise shaping DAC
%
% We can still gain more resolution by shaping the noise psd, i.e. by
% pushing the psd of the quantization noise towards frequencies far above
% Fs/2 (up to 4*Fs/2), while keeping the psd of the signal untouched. We
% will now show that this makes it possible to produce 4-bit quantization
% at 4*Fs with the same resolution as 6-bit quantization at Fs.
%
% We start by quantizing to 6 bits and oversampling by 4.

signal_std=1;
peak=sqrt(signal_std*2);
signal=peak*sin(2*pi*200*(0:1/8000:2));

%quantization
signal_quantized=uquantize(signal,6,2);
snr_quantized=snr(signal,signal_quantized,0)

% Oversampling by 4
signal_pulses=zeros(1,4*length(signal_quantized));
signal_pulses(1:4:4*length(signal_quantized))=...
             signal_quantized;

lp_fir=firpm(300,[0 .24 .26 1],[1 1 0 0]); 
signal_interpolated=filter(lp_fir,1,4*signal_pulses);
pwelch(signal_interpolated,[],[],[],8);

%%
% Let us now perform 4-bit requantization of the interpolated 
% signal, with noise shaping.

N=4;
quantizer_saturation=2; 
q=2*quantizer_saturation/2^N;

delay_memory=0;
signal_requantized=zeros(size(signal_interpolated));
for i=1:length(signal_interpolated)

    u=signal_interpolated(i)-delay_memory;
    % quantization, including dithering
    signal_requantized(i)=uquantize(u,N,quantizer_saturation);
    delay_memory=signal_requantized(i)-u;

end;

signal_sampled_at_4Fs=peak*sin(2*pi*200*(0:1/32000:2));
plot(signal_sampled_at_4Fs(2000:2400)); hold on;
plot(signal_requantized(2151:2550)); hold off;
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% The psd of the requantized signal shows colored quantization noise with
% most of its energy above Fs/2, as a result of noise shaping. Therefore
% the measured SNR does not reflect the actual SNR.

pwelch(signal_requantized,[],[],[],8);
snr_requantized=snr(signal_sampled_at_4Fs,...
                            signal_requantized,300)

%%
% Again, we can SIMULATE the analog part of D/A conversion.

signal_pulses=zeros(1,3*length(signal_requantized));
signal_pulses(1:3:3*length(signal_requantized))=...
                                                        signal_requantized;
hold_impresp=ones(1,3);
signal_staircase=conv(signal_pulses,hold_impresp);

[order,Wn]=cheb1ord(2*pi*3400,2*pi*28600,0.1,60,'s'); 
[Num_LP,Den_LP]=cheby1(order,0.1,Wn,'s'); % relaxed filter

analog_output=filters(Num_LP,Den_LP,signal_staircase,96000);

plot(analog_output(2000:3300));
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% As a result of the analog filtering, the apparent SNR is now only 3 dB
% lower than for our initial 6-bit quantization, i.e. only 3dB less than
% 4-bit quantization at 8 kHz. Again, the lost dBs come from dithering,
% and from the fact that the LP filter specifications are loose.

pwelch(analog_output,[],[],[],24);
signal_sampled_at_12Fs=peak*sin(2*pi*200*(0:1/96000:2));
snr_analog=snr(signal_sampled_at_12Fs,analog_output,4000)
soundsc(analog_output,96000);

%%
% This technique was used in early CD players, when only 14-bit hardware
% D/A converters were available at low cost. By combining oversampling and
% noise shaping (in the digital domain), 14-bit D/A was made comparable to
% 16-bit conventional D/A.

%% 6. Delta-sigma DAC
% The principle of delta-sigma modulation is to apply noise shaping and
% oversampling in such a way that the signal ends up being quantized on 1
% bit. This considerably alleviates the task of the hardware D/A converter
% (which reduces to a switch), and puts most of the load in the digital
% domain. 
% 
% In this proof-of-concept, we will perform 6-to-1 bit requantization,
% using an oversampling ratio of 32 (2^5) and a first-order delta-sigma
% modulators. Real delta-sigma DACs use higher-order modulators, or
% cascades of first-order modulators, and thus do not have to implement
% oversampling of ratio 2^15!
%
% We start by quantizing to 6 bits and oversampling by 32.

signal_std=1;
peak=sqrt(signal_std*2);
signal=peak*sin(2*pi*200*(0:1/8000:1/8)); 

%quantization
N=6;
quantizer_saturation=3;
signal_quantized=uquantize(signal,N,quantizer_saturation);
snr_quantized=snr(signal,signal_quantized,0)

% Oversampling by 32
signal_pulses=zeros(1,32*length(signal_quantized));
signal_pulses(1:32:32*length(signal_quantized))=...
                                                       signal_quantized*32;
lp_fir=firpm(300,[0 1/32-1/100 1/32+1/100 1],[1 1 0 0]);
signal_interpolated=filter(lp_fir,1,signal_pulses);

pwelch(signal_interpolated,[],[],[],64);

%%
% Let us now perform 1-bit requantization with noise shaping. The resulting
% quantized signal is purely binary. It is sometimes referred to as "pulse
% density modulated" (PDM), as the density of its binary transitions is a
% function of the amplitude of the original signal.

N=1;
quantizer_saturation=3; 
q=2*quantizer_saturation/2^N;

signal_requantized=zeros(size(signal_interpolated));
delay_memory=0;
for i=1:length(signal_interpolated)

    u=signal_interpolated(i)-delay_memory;
    % quantization, including dithering
    signal_requantized(i)=uquantize(u,N,quantizer_saturation);
    % saving the internal variable
    delay_memory=signal_requantized(i)-u;

end;

signal_sampled_at_32Fs=peak*sin(2*pi*200*(0:1/256000:1/8));
plot(signal_requantized(2151:3450),'g'); hold on;
plot(signal_sampled_at_32Fs(2000:3300),'linewidth',2); hold off;
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% The psd of the requantized signal shows coloured quantization noise with
% most of its energy above Fs/2, as a result of noise shaping. Therefore
% the measured SNR does not reflect the actual SNR.

pwelch(signal_requantized,[],[],[],64);
snr_requantized=snr(signal_sampled_at_32Fs,signal_requantized,300)

%%
% Again, we can SIMULATE the analog part of D/A conversion. We use an
% "analog sampling frequency" equal to the one we have reached after 32
% times interpolation.

[order,Wn]=cheb1ord(2*pi*3400,2*pi*28600,0.1,60,'s'); 
[Num_LP,Den_LP]=cheby1(order,0.1,Wn,'s'); % relaxed filter
zplane(Num_LP,Den_LP);

analog_output=filters(Num_LP,Den_LP,signal_requantized,256000);

plot(analog_output(2000:3300));
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% As a result of the analog filtering, the apparent SNR is now very close
% to that of our initial 6-bit quantization! This is confirmed by
% listening.

pwelch(analog_output,[],[],[],64);
snr_analog=snr(signal_sampled_at_32Fs,analog_output,4000)
soundsc(analog_output,256000);

%%
% This technique is used in most CD players today, for requantizing 16-bits
% samples to one bit, hence with much higher oversampling ratios.
% Interpolation is therefore performed in several steps, for keeping the
% interpolation filters as simple as possible. Delta-sigma is also the
% heart of the DSD (direct stream digital) coding used in Super-Audio CDs
% (SACD), in which a 1-bit stream is created by the ADC, stored on the CD,
% and directly converted to sound by the DAC. 
% 
% Notice that all the filters used in this proof-of-concept use floating
% point arithmetics. Using fixed-point arithmetics is more complex.

%% Appendix 1: Understanding pwelch
% Since it is used al over this script, it is important to understand what
% MATLAB's |pwelch| function shows, and why. 
%
% |pwlech| plots an estimate of the power spectral density S(f) of a signal
% x(t) by computing the power spectrum P(f) of its samples x(n). The
% estimate of P(f) is obtained by computing the periodogram (also called
% the Welch estimator) of x(n).
%
% The input signal is divided into frames (8 frames by default), and for
% each frame a power spectrum estimate (a periodogram) is computed as
% ||fft||²/N, where N is the number of samples in the frame. The average
% periodogram computes the ensemble average of these individual
% periodograms. 
%
% The PSD S(f) is the Fourier Transform of the autocorrelation function
% phi(t) of the signal x(t), while the Power Spectrum P(f) is the (discrete
% time) Fourier Transform of the autocorrelation sequence phi(n) of x(n).
% Since the values of phi(n) are samples of phi(t), and assuming the
% sampling frequency Fs was chosen so as to avoid aliasing, then P(f) is
% simply equal to S(f)*Fs in the range [-Fs/2,Fs/2]. As a result, |pwelch|
% returns (and shows) P(f)/Fs, i.e. ||fft||²/N/Fs.
%
% What is more, when used for computing the PSD of real signal, |pwelch|
% returns by default the onesided PSD, i.e. the PSD over [0,Fs/2], which is
% itslef normalized so that its integral over [0,Fs/2] yields the expected
% power of the signal. In other words, |pwelch| computes (and shows)
% 2*||fft||²/N/Fs.
%
% Consequently, claiming Fs=2 as an input argument for |pwelch| shows
% ||fft||²/N in the range [0,1], as did the good old, but now obsolete,
% |psd| function. The integral of this function over [0,1] (i.e., the
% power, or variance, of the signal) is thus the average of the psd shown
% on the plot. 
% Notice, also, that estimating this average by eye is not easy, since
% psds are plotted in dB. When the signal is white noise, however, its
% variance can be read directly as the average of the psd.
%
% As an example, let us plot the psd of Gaussian white noise with zero mean
% and variance=100, and claim Fs=2. |pwelch| shows an average value of 20
% dB.

signal=randn(1,10000)*10;
pwelch(signal,[],[],[],2);

%%
% (in reality, |pwelch| also uses weighting windows and compensates for
% their effect on the psd estimate)

%% 
% Claiming Fs=10 divides the psd by 10, i.e, 10 dB.
pwelch(signal,[],[],[],20);

%%
% Notice finally that a sinusoid with amplitude A will appear via |pwelch|
% as a spectral line equal to 2 A²(N/2)²/(N Fs), i.e., A²(N/2)/Fs. For a
% sine wave of amplitude 10 whose periodogram is computed on 128 points,
% claiming Fs=2 gives 10*log10(100*128/2/2)=35 dB. In order to get this
% value on the screen, we need to force frames to be weighted by a
% rectangular window (i.e., not weihgted at all).

signal=10*sin(2*pi*100*(0:1/1000:1)); % 100 Hz at Fs= 1000 Hz
pwelch(signal,rectwin(128),[],[],2);

%%
% If we let |pwelch| use its default window (a Hamming window), it will not
% show exactly the expected value, because the compensation factor it
% applies on the psd assumes that this signal is random. As a result, the
% values of spectral lines shown by this function should be only taken as
% an gross indication. MATLAB actually provides another function, which returns
% the true psd of sinusoids: A²/2. In our case : 10*log10(100*/2)=17 dB

msspectrum(spectrum.welch,signal);
