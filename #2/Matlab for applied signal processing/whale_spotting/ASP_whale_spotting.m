%% Chapter 6 - How can marine biologists track sperm whales in the oceans?
% This is a companion file to the book "Applied Signal Processing", 
% by T.Dutoit and F. Marques, Springer 2008.
% 
% It is supposed to be run cell-by-cell, using the cell mode of 
% MATLAB 6 and later versions. Search for "what are cells" in the 
% product help of MATLAB to know how to use them.
% 
% This file uses the SIGNAL_PROCESSING toolbox of MATLAB.

%%
% In this script, we analyze sperm whale clicks received by several
% hydrophones (Section 1), and use them to find the position of the
% cetacean. We start, in Section 2, by processing the clicks with the
% (non-linear) Teager-Kaiser operator, so as to increase their
% signal-to-noise ratio. We then test cross-correlation estimation
% algorithms, in Section 3, to compute the time difference of arrival
% (TDOA) between pairs of hydrophones, and compare them in Section 4, to an
% efficient adaptive estimation algorithm based on least means squares
% minimization. We conclude in Section 5 by a simple multilateration
% algorithm, which uses two TDOAs to estimate the position of the sperm
% whale.
%
% Copyright V. Kandia, Y. Stylianou, T. Dutoit (2007)

clear all;
set(0,'defaultFigureColor','w');

%% 1. Sperm whale sounds
% Throughout this script, we will work on hydrophone signals collected at
% the Atlantic Undersea Test and Evaluation Center (AUTEC), Andros Island,
% Bahamas, and made available by the Naval Undersea War-fare Center (NUWC).
% 
% These signals were provided to the 2nd International Workshop on
% Detection and Localization of Marine Mammals using Passive Acoustics,
% held in Monaco November 16–18, 2005 and can be obtained from their website.
% The Andros Island area has over five hundred square nautical miles of
% ocean that are simultaneously monitored via 68 broad-band hydrophones.
% The distance between most hydrophones is 5 nautical miles (about 7.5 km).
% We will work with the sounds received by hydrophones I, G, and H (which
% we will rename as #1, #2, and #3 here), for about 30 seconds, recorded
% with a DAT at Fs=48 kHz. 
%
% Plotting the signals received by all three hydrophones shows variable
% attenuations and delays between clicks, as well as variable noise levels.
% Each direct path click is followed by an echo due to surface reflection. 
% Obviously, clicks arrive at hydrophone #2 much earlier than at the other
% hydrophones. The amplitude of the signal on hydrophone #2 is also higher,
% and the noise level is lower. Moreover, the echo in signal #2 is much
% more delayed than in signals #1 and #3. All this suggests that the sperm
% whale is in the vicinity of hydrophone #2.

[hydrophone1,Fs]=wavread('hydrophone1.wav');
[hydrophone2,Fs]=wavread('hydrophone2.wav');
[hydrophone3,Fs]=wavread('hydrophone3.wav');

clear ax;
time=(0:length(hydrophone1)-1)/Fs;
ax(1)=subplot(3,1,1); plot(time,hydrophone1)
ylabel('Amplitude');
title('hydrophone #1');
hold on;
ax(2)=subplot(3,1,2); plot(time,hydrophone2)
ylabel('Amplitude');
title('hydrophone #2');
ax(3)=subplot(3,1,3); plot(time,hydrophone3)
xlabel('Time (s)'); ylabel('Amplitude');
title('hydrophone #3');
linkaxes(ax,'x');

%%
% Let us listen to the first two clicks at hydrophone #1. 

soundsc(hydrophone1(7.5*Fs:10.5*Fs),Fs);

%% 
% Since the delay between clicks is larger than the interclick interval, we
% will first apply a fixed time-shift between signals, so as to
% approximately align corresponding clicks. We shift the first signal by
% 7.72 ms to the left, the second signal by 5.4 s to the left, third signal
% by 7.58 s to the left, and keep 25 s of each signal. We will later
% account for these artificial shifts. 

n_samples=25*Fs;
shift_1=fix(7.72*Fs);
shift_2=fix(5.40*Fs);
shift_3=fix(7.58*Fs);
[hydrophone1,Fs]=wavread('hydrophone1.wav',[1+shift_1 shift_1+n_samples]);
[hydrophone2,Fs]=wavread('hydrophone2.wav',[1+shift_2 shift_2+n_samples]);
[hydrophone3,Fs]=wavread('hydrophone3.wav',[1+shift_3 shift_3+n_samples]);

clear ax;
time=(0:length(hydrophone1)-1)/Fs;
ax(1)=subplot(3,1,1); plot(time,hydrophone1)
ylabel('Amplitude');
title('hydrophone #1');
hold on;
ax(2)=subplot(3,1,2); plot(time,hydrophone2)
ylabel('Amplitude');
title('hydrophone #2');
ax(3)=subplot(3,1,3); plot(time,hydrophone3)
xlabel('Time (s)'); ylabel('Amplitude');
title('hydrophone #3');
linkaxes(ax,'x');

%%
% Zooming on a single click shows that
% whale clicks are not as simple as Dirac samples, which makes the visual
% estimation of TDOAs not very precise, and their automatic estimation not
% trivial.

set(gca,'xlim',[6.66 6.73]);

%%
% Let us compute the PSD of the background noise in hydrophone #3, for
% later use. It appears that the noise is pink, with 40 dB more power
% around f=0 Hz than around f=24000 Hz.

real_noise=hydrophone3(1.26e4:1.46e4);
real_noise_std=std(real_noise)

clf;
pwelch(real_noise,[],[],[],Fs);

%%
% Notice that the data used here does not reveal LF interference, which may
% however occur in measurements and make the estimation of TDOA
% still more complex.

%% 2. Teager-Kaiser filtering
% In order to better understand the Teager-Kaiser non-linear filter, we
% first apply it to synthetic signals, composed of clicks (assimilated to
% Dirac impulses), background noise (white and pink noise), and an
% interference component (modeled as a sinusoid). We then show that, 
% even though the effect of a general non-linear filter on a sum of signals 
% is not the sum of its outputs to isolated input signals, the effect of 
% the TK operator on synthetic clicks with background noise and sinusoidal 
% interference can still somehow be analyzed in terms of its effect on 
% isolated components. We finally apply the TK operator on real sperm whale 
% sounds, and compare it to a simple linear high-pass filter.

%% 2.1 Dirac pulse input
% We start with the response of the TK operator to a Dirac impulse at
% n=4000, with the same order of magnitude as that of the clicks we found
% in the previous Section. The output is the input impulse, squared.  

clf;
click = [zeros(1,3999) 0.1 zeros(1,4000)];

% Applying TK filter
L = length(click);
click_response = click(2:L-1).^2-click(1:L-2).*click(3:L);
click_response = [click_response(1) click_response ...
    click_response(L-2)]; % lenght(output)=length(input)

subplot(2,1,1)
plot(click);
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(2,1,2)
plot(click_response);
xlabel('Time (samples)'); ylabel('Amplitude');

%% 2.2 Sinusoidal input
% Sinusoidal interference signals of the form A cos(n phi0) produce a constant
% output A²sin²(phi0). 
% To check this we generate 8000 samples of a chirp cos(n phi_0(n)) with
% phi_0(n)=phi_max n/8000 and phi_max=0.1 which corresponds to a
% frequencies from 0 to Fs phi_max/(2 pi) = 763 Hz.   
%
% *MATLAB function involved:*
% 
% * |y = teager_kaiser(x)| applies the Teager-Kaiser energy operator to
% input signal |x|.

Fs=48000;
phi_max=0.1;
test=chirp(0:1/Fs:7999/Fs,0,7999/Fs,Fs*phi_max/2/pi); 
test_response=teager_kaiser(test);

clf;
subplot(2,1,1)
plot(test);
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(2,1,2)
plot(test_response);
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% The value of the output is close to zero when phi0
% is small, which is be the case for a typical LF interference component in
% sperm whale click measurements.

interference = 0.01*cos (0.01*(0:7999)); % A²sin²(0.01)~1e-8
interference_response=teager_kaiser(interference);

subplot(2,1,1)
plot(interference);
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(2,1,2)
plot(interference_response);
set(gca,'ylim',[0 1.2e-8]);
xlabel('Time (samples)'); ylabel('Amplitude');

%% 2.3 White noise input
% TK filtering of Gaussian white noise N(0,sigma) produces noise with
% mean of the order of sigma² and standard deviation proportional to
% sigma². We check this on white noise with linearly increasing
% standard deviation. 

test = (1/8000:1/8000:1).*randn(1, 8000);
test_response = teager_kaiser(test);

subplot(2,1,1);
plot(test);
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(2,1,2);
plot(test_response);
set(gca,'ylim',[-3 5]);
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% The output noise in still white.

noise = randn(1, 8000);
noise_response = teager_kaiser(noise);

close(1);
subplot(2,1,1)
pwelch(noise,[],[],[],Fs);
subplot(2,1,2)
pwelch(noise_response,[],[],[],Fs);

%% 
% But it is not Gaussian.

clf; clear ax;
ax(1)=subplot(2,1,1);
hist(noise,30);
ax(2)=subplot(2,1,2);
hist(noise_response,60);
set(gca,'xlim',[-10 10]);

%% 2.4 Pink noise input
% Let us now test the TK filter on pink noise similar to the one we found
% on real hydrophone signals in Section 1. We first create a filter with
% linear frequency response from -30 dB at f=0 to -65 dB at f=Fs/2. 

clf;
w=0:0.1:1;
a_dB=-35*w-30; 
a=10.^(a_dB/20);
[B,A] = fir2(20,w,a);
freqz(B,A,512,48000);

%%
% Applying this filter to white noise N(0,1) produces a realistic pink 
% noise component. 

pink_noise=filter(B,A,noise); 
pwelch(pink_noise,[],[],[],Fs);
pink_noise_std=std(pink_noise)

%%
% TK filtering of pink noise produces pink noise, whose standard deviation
% is still proportional to the square of that of the noise (although with
% a smaller the proportionality factor than for white input noise).  

test=100*(1/8000:1/8000:1).*pink_noise;
test_response = teager_kaiser(test);

subplot(2,1,1);
plot(test);
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(2,1,2);
plot(test_response);
set(gca,'ylim',[-3 5]);
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% The output noise is whiter; the PSD of the input noise has been decreased
% by about 50 dB in low frequencies and by about 30 dB in high frequencies.

pink_noise_response = teager_kaiser(pink_noise);

close(1);
subplot(2,1,1)
pwelch(pink_noise,[],[],[],Fs);
subplot(2,1,2)
pwelch(pink_noise_response,[],[],[],Fs);

%% 2.5 Complex input
% We now apply the TK operator to the complete synthetic signal, obtained
% by summing all three components: click, interference, and pink noise. 

input=click+interference+pink_noise;
output_tk=teager_kaiser(input);

%%
% Thanks to the TK operator, the impulse is highlighted in the signal. The LF
% signal is removed, and the pink noise is very much attenuated relatively
% to the impulse. The SNR is therefore highly increased. 

clear ax;
ax(1)=subplot(2,1,1);
plot(input);
xlabel('Time (samples)'); ylabel('Amplitude');
ax(2)=subplot(2,1,2);
plot(output_tk);
xlabel('Time (samples)'); ylabel('Amplitude');
linkaxes(ax,'x');

%%
% Notice that the TK operator did not smear the impulsive part of the input
% waveform. 

set(gca,'xlim',[3980 4020]);

%%
% Notice also that the effect of the TK filter on the SNR depends on the
% original value of the SNR. As a matter of fact, TK squares the amplitude
% A of the impulse and the standard deviation sigma of the noise. 
%
% Let us characterize the "visibility" of the pulse by the ratio
% A/(3*sigma), since noise samples take most of their values in [-3*sigma,
% 3*sigma]. 
% If this ratio is significantly higher than 1 (in the previous plot it was
% 0.1/0.033=3),  the visibility of the pulse in the TK filtered signal,
% A²/(9*sigma²), is higher than in the original signal. If A/(3 sigma)
% approaches 1 but is still higher, the impulse may not be visible in the
% original signal, but more easily detected in the output signal.   
%
% To check this, we multiply the standard deviation of the noise by two. 

input=click+2*pink_noise+interference;
output_tk=teager_kaiser(input);

clear ax;
ax(1)=subplot(2,1,1);
plot(input);
xlabel('Time (samples)'); ylabel('Amplitude');
ax(2)=subplot(2,1,2);
plot(output_tk);
xlabel('Time (samples)'); ylabel('Amplitude');
linkaxes(ax,'x');

%%
% Summarizing, the TK operator has a high-pass filtering effect on
% background pink noise (and it squares its standard deviation), a notch
% effect on possible low-frequency interference (while a DC component
% appears), and amlost no effect on Dirac pulses (aprat from squaring their
% amplitude): it acts as a signal-dependent filter and amplifier. 
% (In practice, it also produces cross terms when several components are
% summed at the input. We did not examine them here, as they are not
% essential for the kind of signal we deal with.)

%% 2.6 Hydrophone signal
% Finally we apply the TK operator to a real click taken from a hydrophone
% signal, as examined in Section 1, and compare its effect to that of a linear filter.
% We create a linear-phase high-pass symmetric FIR filter whose frequency
% response has approximately the same effect on the pink background noise as the
% TK filter: a frequency-linear attenuation, from -50 dB close to f=0 
% -30 dB close to f=Fs/2. 

clf;
w=0:0.1:1;
a_dB=20*w-50;
a=10.^(a_dB/20);
N=30;
[B,A] = fir2(N,w,a);
freqz(B,A,512,Fs);

%%
% The length of its impulse response is 31 (its order is 30). It therefore 
% has a delay of 15 samples.  

plot(B);
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% Obviously the TK operator does a very good job at highlighting the click.
% While the linear filter has a similar effect on background noise, it also
% decreases the power of the click, thereby increasing the overall SNR. 

input=wavread('hydrophone3.wav',[465000 466650]);

output_tk=teager_kaiser(input);
% Taking the filter dealy into account
output_lin=[filter(B,A,input(16:end)) ; zeros(15,1)];

ax(1)=subplot(3,1,1);
plot(input);
ylabel('Amplitude');
ax(2)=subplot(3,1,2);
plot(output_tk);
set(gca,'ylim',[-0.01 0.01]);
ylabel('Amplitude');
ax(3)=subplot(3,1,3);
plot(output_lin);
set(gca,'ylim',[-0.01 0.01]);
xlabel('Time (samples)'); ylabel('Amplitude');
linkaxes(ax,'x');

%% 3. TDOA estimation using generalized cross-correlation
% We will now show how generalized cross correlation (GCC) can be used for
% estimating TDOAs. We start again by applying it to simple synthetic signals,
% and end with real hydrophone signals.

%% 3.1 Synthetic input
% We first create a simple synthetic signal, as in Section 2 but with more
% severe interference and noise, and a second signal, with other noise
% samples, a phase-shifted interference signal, and an attenuated click,
% delayed by 60 samples. 

click = [zeros(1,499) 0.1 zeros(1,500)];
interference = 0.03*cos (0.01*(0:999));
noise = 2*randn(1, 1000);
w=0:0.1:1;
a_dB=-35*w-30; 
a=10.^(a_dB/20);
[B,A] = fir2(20,w,a);
pink_noise=filter(B,A,noise); 

signal_1=click+interference+pink_noise;

click_2=circshift(click'*0.5,60)';
interference_2 = 0.02*cos (pi/4+0.01*(0:999));
noise_2 = randn(1, 1000); % new noise samples
pink_noise_2=filter(B,A,noise_2);

signal_2=click_2+interference_2+pink_noise_2;

subplot(2,1,1);
plot(signal_1);
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(2,1,2);
plot(signal_2);
xlabel('Time (samples)'); ylabel('Amplitude');

save synthetic_signals signal_1 signal_2 % for use in Section 3

%% 
% Estimating the time difference of arrival (TDOA) between these signals
% with the standard cross-correlation does not provide the expected value
% (-60 samples), as the total cross-correlation function is dominated by
% correlation between the sinusoidal interference components.
% Notice the (biased) cross-correlation is obtained via FFT-IFFT, after
% making sure that the number of FFT samples is high enough for linear
% convolution to be identical to circular convolution (as obtained by
% multiplying FFTs). 
%
% Notice that the linear phase component in the HF part of the cross PSD,
% where the pink noise is weak, reveals a delay. But since the amplitude of
% the noise and of the LF interference strongly dominate that of the HF
% part of the cross PSD, this delay cannot be observed in the cross
% corre-lation function.    

M = min(length(signal_1),length(signal_2));
NFFT = 2*M-1; % so that linear convolution = circular convolution 
x1 = signal_1 - mean(signal_1);
x2 = signal_2 - mean(signal_2);
X1 = fft(x1,NFFT);
X2 = fft(x2,NFFT);
S_x1x2 = X1.*conj(X2);
phi_x1x2 = ifft(S_x1x2);
% re-arranging the IFFT
standard_xcorr = [phi_x1x2(NFFT-M+2:NFFT) phi_x1x2(1:M)];  
[val,ind]=max(standard_xcorr);
TDOA_standard_xcorr = ind-M  % in samples

clf;
subplot(2,1,1);
plot((0:1/M:(M-1)/M), 20*log10(abs(S_x1x2(1:M))) );
set(gca,'ylim',[-50 50]);
xlabel('Frequency (x pi rad/sample)'); ylabel('|S_x_1_x_2| (dB)');
subplot(2,1,2);
plot((0:1/M:(M-1)/M), unwrap(angle(S_x1x2(1:M))) );
xlabel('Frequency (x pi rad/sample)'); ylabel('arg(S_x_1_x_2) (rad)');

%%
clf;
plot((-M+1:M-1),standard_xcorr); grid;
xlabel('Time lag (sample)'); ylabel('Amplitude');

%%
% Using the phase transform version of the generalized cross-correlation
% removes the LF interference and gives the same importance to all
% frequency bands in the phase spectrum.  It produces a more prominent
% maximum, which leads to a correct estimate of the TDOA 

PT= S_x1x2 ./ max(abs(S_x1x2),eps);
phi_x1x2 = ifft(PT);
% re-arranging the IFFT
phase_transform = [phi_x1x2(NFFT-M+2:NFFT) phi_x1x2(1:M)]; 

[val,ind]=max(phase_transform);
TDOA_phase_transform = ind-M

%%
clf;
subplot(2,1,1);
plot((0:1/M:(M-1)/M), 20*log10(abs(PT(1:M))) );
set(gca,'ylim',[-50 50]);

xlabel('frequency (x pi rad/sample)'); ylabel('|PT| (dB)');
subplot(2,1,2);
plot((0:1/M:(M-1)/M), unwrap(angle(PT(1:M))) );
xlabel('frequency (x pi rad/sample)'); ylabel('arg(PT) (rad)');


%%
clf;
plot((-M+1:M-1),phase_transform); grid;
xlabel('Time lag (sample)'); ylabel('Amplitude');


%%
% Applying the TK operator obviously increases the SNR.

signal_1_tk=teager_kaiser(signal_1);
signal_2_tk=teager_kaiser(signal_2);

subplot(2,1,1);
plot(signal_1_tk);
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(2,1,2);
plot(signal_2_tk);
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% As a result of the suppression of the sinusoidal component, GCC produces
% an accurate result. 
%
% *MATLAB function involved:*
% 
% * |gcc(z1, z2, flag)| computes the generalized cross correlation (GCC)
% between signals z1 and z2, from FFT/IFFT, as specified in (Knapp & Carter
% 1976). |[flag]| makes it possible to choose the type of cross-correlation:
%  Standard cross correlation if |flag='cc'|; Phase transform:if |flag='phat'|

standard_xcorr_tk = gcc(signal_1_tk, signal_2_tk, 'cc');
[val,ind]=max(standard_xcorr_tk);
TDOA_standard_xcorr_tk = ind-M  % in samples

clf;
plot((-M+1:M-1),standard_xcorr_tk); grid;
xlabel('Time lag (sample)'); ylabel('Amplitude');

%%
% Using the phase transform on TK-filtered data still produces a more
% prominent maximum. 

phase_transform_tk = gcc(teager_kaiser(signal_1), teager_kaiser(signal_2),  'phat');
[val,ind]=max(phase_transform_tk);
TDOA_phase_transform_tk = ind-M

clf;
plot((-M+1:M-1),phase_transform_tk); grid;
xlabel('Time lag (sample)'); ylabel('Amplitude');

%% 3.2 Hydrophone signals
% We now apply GCC to real hydrophone signals (notice we . Clicks are separated by
% about 400 samples.

Fs=48000;
shift_2=fix(5.40*Fs);
shift_3=fix(7.58*Fs);

signal_1=wavread('hydrophone2.wav',[460000+shift_2 ...
    460000+shift_2+2000]);
signal_2=wavread('hydrophone3.wav',[460000+shift_3 ...
    460000+shift_3+2000]);

subplot(2,1,1);
plot(signal_1);
ylabel('Amplitude');
subplot(2,1,2);
plot(signal_2);
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% Again, applying the TK operator has a positive effect on the SNR.

signal_1_tk=teager_kaiser(signal_1);
signal_2_tk=teager_kaiser(signal_2);

subplot(2,1,1);
plot(signal_1_tk);
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(2,1,2);
plot(signal_2_tk);
xlabel('Time (samples)'); ylabel('Amplitude');

%% 
% Applying the Phase Transform to these signals delivers an estimate of the
% TDOA: 404 samples (negative, since the click hits hydrophone #1 before
% hydrophone #2). 

phase_transform_tk  = gcc(signal_1_tk, signal_2_tk,  'phat');
[val,ind]=max(phase_transform_tk );

M = min(length(signal_1),length(signal_2));
TDOA_phase_transform_tk= ind-M

clf;
plot((-M+1:M-1),phase_transform_tk); grid;

%% 4 TDOA estimation using least-mean squares
% In this section we will apply the Least Mean Square (LMS) adaptive
% approach to TDOA estimation, again first on synthetic signals (with
% imposed TDOA), and then on real hydrophone signals.  

%% 4.1 Synthetic input
% We start with the same synthetic signals as in Section 2, and apply
% TK preprocessing. 

load synthetic_signals signal_1 signal_2 % computed in Section 2

signal_1_tk=teager_kaiser(signal_1);
signal_2_tk=teager_kaiser(signal_2);

subplot(2,1,1);
plot(signal_1_tk);
xlabel('Time (samples)'); ylabel('Amplitude');
subplot(2,1,2);
plot(signal_2_tk);
xlabel('Time (samples)'); ylabel('Amplitude');

%%
% We run the adaptive filtering algorithm for TDOA estimation within [-600,600]
% (samples) and step |mu|=0.01 (for full scale signals in [-1,+1]).

% Normalizing max signal amplitudes to +1
x1=signal_1_tk/max(signal_1_tk);
x2=signal_2_tk/max(signal_2_tk);

% LMS initialization
M = 600; % max value of the estimated TDOA
x1c = zeros(M,1);
x2c = zeros(M,1);
u = zeros(2*M,1);
u(M/2) = 1;
N = length(x1);
e = zeros(1,N);
tdoa = zeros(1,N);
peak = zeros(1,N);
mu = 0.01; % LMS step

% LMS loop
for n=1:N
    
    x1c = [x1(n);x1c(1:length(x1c)-1)];
    x2c = [x2(n);x2c(1:length(x2c)-1)];
    x = [x1c;x2c];
    
    e(n) = u'*x;
    u = u-mu*e(n)*x;
    u(M/2) = 1; %forcing g2 to an impulse response at M/2
    u = u/norm(u); %forcing ||u|| to 1
    
    [peak(n),ind] = min(u(M+1:end));
    peak(n)=-peak(n); % find the value of the (positive) impulse
    TDOA(n) = ind-M/2;
    
end

% Estimated TDOA as a function of time, with values of the peak in h1
subplot(2,1,1);
plot(TDOA);
xlabel('Time (samples)'); ylabel('TDOA (samples)');
subplot(2,1,2);
plot(peak);
xlabel('Time (samples)'); ylabel('peak');

%%
% It appears that the TDOA estimate is wrong at the beginning of the frame,
% as the adaptive filter must have started processing a click to start
% converging. The value of the peak in h1 (the impulse response between the
% source and the first signal) is very low for these wrong TDOA values. The
% best TDOA estimate is the one that produced the most prominent peak in
% h1. This estimate is correct: -60 samples.
% After the clicks have been processed, the peak in h1 starts decreasing again, but
% the value of TDOA remains correct. 

[val,ind]=max(peak);
Best_estimate_TDOA=TDOA(ind)

%% 4.2 Hydrophone signals
% Applying the same algorithm to hydrophone signals #1 and #2 is
% straightforward (notice we first time-shift the signals, as in Section 1).

Fs=48000;
n_samples=25*Fs;
shift_1=fix(7.72*Fs);
shift_2=fix(5.40*Fs);
[hydrophone1,Fs]=wavread('hydrophone1.wav',[1+shift_1 shift_1+n_samples]);
[hydrophone2,Fs]=wavread('hydrophone2.wav',[1+shift_2 shift_2+n_samples]);

% Subsampling by 6 (-> Fs'=8 kHz), for decreasing computational cost
signal_1=resample(hydrophone1,1,6);
signal_2=resample(hydrophone2,1,6);

% Applying Teager-Kaiser filter
signal_1_tk=teager_kaiser(signal_1);
signal_2_tk=teager_kaiser(signal_2);

clear ax;
ax(1)=subplot(2,1,1);
plot((0:length(signal_1_tk)-1)/(Fs/6),signal_1_tk);
ylabel('Amplitude');
ax(2)=subplot(2,1,2);
plot((0:length(signal_2_tk)-1)/(Fs/6),signal_2_tk);
xlabel('Time (s)'); ylabel('Amplitude');
linkaxes(ax,'x');

%%
% The estimated TDOA change by a few tens of milliseconds in our 25-seconds
% recording. This shows that the sperm whale moved during the recording. 
% Again, the first TDOA estimates are not significant. 
% Plotting our TODAs with the ones obtained by visual inspection of the
% signals shows that the LMS estimation provided accurate results.
%
% *MATLAB function involved:*
% 
% * |[tdoa, peak] = TDOA_LMS(x1, x2, max_tdoa, mu)| is an implementation of 
% (Benesty 2000) for source localization. It returns the time difference of
% arrival, in samples, between signals |x1| and |x2|, as a function of
% time, together with the value of the peak found in the estimate of the
% main propagation path for signal |x1|. This value can be used as a
% confidence value for the TDOA. |max_tdoa| is the maximum TDOA value
% returned and |mu| is the step weight.

TDOA_12 = TDOA_LMS(signal_1_tk,signal_2_tk,600,0.01);
TDOA_12 = TDOA_12*6; % taking downsampling into account
% Accounting for preliminary time-shifts
TDOA_12 = (TDOA_12 + shift_1-shift_2)/Fs; 

load clicks; % clicks_1, clicks_2, and clicks_3 give the positions 
             % of clicks in the respective hydrophone signals, as estimated
             % by visual inspection. 
TDOA_12_visual=(clicks_1-clicks_2)/Fs;

clf;
plot((0:length(TDOA_12)-1)/(Fs/6),TDOA_12);
hold on;
stairs((clicks_2-shift_2)/Fs,TDOA_12_visual,'--r'); 
xlabel('Time (s)'); ylabel('TDOA (s)');
legend('Estimated TDOA_1_2','Visual TDOA_1_2')
set(gca,'ylim',[2.324 2.33]);
hold off;

%%
% We repeat the operation for hydrophone signals #2 and #3.

shift_3=fix(7.58*Fs);
[hydrophone3,Fs]=wavread('hydrophone3.wav',[1+shift_3 shift_3+n_samples]);

signal_3=resample(hydrophone3,1,6);
signal_3_tk=teager_kaiser(signal_3);

[TDOA_23,peak] = TDOA_LMS(signal_2_tk,signal_3_tk,600,0.01);
TDOA_23 = TDOA_23*6; % taking downsampling into account
% Accounting for preliminary time-shifts
TDOA_23 = (TDOA_23 + shift_2-shift_3)/Fs;

TDOA_23_visual=(clicks_2-clicks_3)/Fs;

clf;
plot((0:length(TDOA_23)-1)/(Fs/6),TDOA_23);
hold on;
stairs((clicks_2-shift_2)/Fs,TDOA_23_visual,'--r'); 
xlabel('Time (s)'); ylabel('TDOA (s)');
legend('Estimated TDOA_2_3','Visual TDOA_2_3')
set(gca,'ylim',[-2.2 -2.18 ]);
hold off;

save tdoas TDOA_12 TDOA_23 % for use in Section 5

%% 5. Multilateration
% In this final Section, we use the TDOAs estimated from real hydrophones
% signals in Section 4, for estimating the position of the sperm whale on a
% 2D map. We implement a rudimentary multilateration system, on the basis
% of 2 TDOAs only. 
% We first show the position of the hydrophones on the map, from their
% (x,y) coordinates, and then plot the hyperbolas corresponding to each
% pair of TDOA, assuming the speed of sound in water is 1510 m/s. Notice
% we drop the first TDOA estimates, which the LMS algorithm used for
% converging, and only use 10 estimates for our plot. 
%
% *MATLAB function involved:*
% 
% * |plot_hyp(x1,y1,x2,y2,a,color)| plots a hyperbola whose foci 
% are (x1,y1) and (x2,y2), and whose semi-major axis is |a|

hydrophone_pos = [14318.86    -16189.18    % hydrophone 1: x,y [meters]
                  10658.04    -14953.63    % hydrophone 2: x,y
                  12788.99    -11897.12];  % hydrophone 3: x,y
plot(hydrophone_pos(:,1),hydrophone_pos(:,2),'o');
text(hydrophone_pos(1,1),hydrophone_pos(1,2),'  1');
text(hydrophone_pos(2,1),hydrophone_pos(2,2),'  2');
text(hydrophone_pos(3,1),hydrophone_pos(3,2),'  3');
xlabel('x (m)'); ylabel('y (m)');
hold on;

load 'tdoas' TDOA_12 TDOA_23

% Keeping 10 TDOA estimates 
TDOA_12=TDOA_12(6000:length(TDOA_12)/10:end);
TDOA_23=TDOA_23(6000:length(TDOA_23)/10:end); 

sound_speed=1510; % m/s
for i=1:length(TDOA_12)
   plot_hyp(hydrophone_pos(1,1), hydrophone_pos(1,2), ...
            hydrophone_pos(2,1), hydrophone_pos(2,2), ...
            -TDOA_12(i)*sound_speed/2,'b');
   plot_hyp(hydrophone_pos(2,1), hydrophone_pos(2,2), ...
            hydrophone_pos(3,1), hydrophone_pos(3,2), ...
            -TDOA_23(i)*sound_speed/2,'r');
end;
hold off;

%%
% Since the maximum speed of a whale is close to 30 km/h (i.e. 8 m/s), it
% obviously could not move much around hydrophone #2 in the 25 s of data we
% have used for this proof-of-concept. 
% Zooming around hydrophone #2 reveals its trajectory : it moved from
% (10837.5,-14822) to (10844,-14847), i.e. by about 26 m in 25 seconds,
% i.e. at a speed of about 1 m/s. 

clf;
axis([10835, 10847, -14856, -14816]);
grid;
xlabel('x (m)'); ylabel('y (m)');
hold on;
for i=1:length(TDOA_12)-1
    plot_hyp(hydrophone_pos(1,1), hydrophone_pos(1,2), ...
             hydrophone_pos(2,1), hydrophone_pos(2,2), ...
             -TDOA_12(i)*sound_speed/2,'b');
    plot_hyp(hydrophone_pos(2,1), hydrophone_pos(2,2), ...
             hydrophone_pos(3,1), hydrophone_pos(3,2), ...
             -TDOA_23(i)*sound_speed/2,'r');
pause(0.5);
end;
hold off;

distance = sqrt((10844-10837.5)^2+(14847-14822)^2)
