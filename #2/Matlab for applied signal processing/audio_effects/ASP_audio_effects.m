%% Chapter 5 - How does an audio effects processor perform pitch-shifting? 
% This is a companion file to the book "Applied Signal Processing", 
% by T.Dutoit and F. Marques, Springer 2008.
% 
% It is supposed to be run cell-by-cell, using the cell mode of 
% MATLAB 6 and later versions. Search for "what are cells" in the 
% product help of MATLAB to know how to use them.
% 
% This file uses the SIGNAL_PROCESSING toolbox of MATLAB.
 
%%
% In this script, we show that STFT-based signal processing provides very
% flexible tools for audio signal processing. We start be radically
% modifying the phases of the STFT, to create a robotization effect
% (Section 1). Then we examine the MATLAB implementation of a phase vocoder,
% with and without vertical phase locking (Section 2). We conclude the Chapter
% with two pitch-scale modification algorithms, which use the phase vocoder
% in various aspects (Section 3).
%
% Copyright T. Dutoit, J. Laroche (2007)

set(0,'defaultFigureColor','w');

%% 1. STFT-based audio signal processing 
% We first implement a basic STFT-based processing scheme, using the same
% weighting window for both analysis and synthesis, and with very simple
% processing of the intermediate DFTs. 

%% 1.1 Weighting windows
% Since weighting windows play an important part in STFT-based signal
% processing, it is interesting check their features first. MATLAB proposes
% a handy tool for that: the |wvtool| function. Let us use it for comparing
% the rectangular, Hanning, and Blackman windows. Clearly, the spectral
% leakage of the Hanning and Blackman windows are lower than those of the
% rectangular and sqrt(Hanning) windows. By zooming on the the main
% spectral lobes, one can check that this is compensated by a higher lobe
% width : the spectral width of the rectangular window is half that of the
% Hanning window, and a third of that of the Blackman window.

N=100;
wvtool(boxcar(N),hanning(N),blackman(N));

%% 1.2 The Constant OverLap-Add (COLA) constraint
% Let us now examine the operations involved in weighted overlap-add
% (WOLA). When no modification of the analysis frames is performed, one
% can to the least expect the output signal to be very close to the input
% signal. In particular, a constant input signal should result in a
% constant output signal. This condition is termed as Constant OverLap-Add
% (COLA) constraint. It is strictly met for specific choices of the
% weighting window and of the window shift, as we shall see. We now check
% this on a chirp signal.

Fs=8000;
input_signal=chirp((0:Fs)/Fs,0,1,2000)';
soundsc(input_signal, Fs);
specgram(input_signal,1024,Fs,256);

%%
% We actually use the periodic version of the Hanning window, which better
% meets the COLA constraint for various values of the frame shift than the
% default symmetric version. It is easy to plot the SNR of the
% analysis-synthesis process, for a variable frame shift between 1 and 512
% samples (and a frame length of 512 samples). 
% The resulting SNR is indeed very high for values of the frame shift equal
% to N/4, N/8, N/16, ..., 1. These values are said to meet the COLA
% contraint for the squared Hanning window. In practice, values of the
% frame shift that do not meet the COLA constraint can still be used,
% provided the related SNR is higher than the local signal-to-mask ratio
% (SMR) due to psycho-acoustic effects (see Chapter 3). Using the highest
% possible value for the frame shift while still reaching a very high SNR
% minimizes the computational load of the system.
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

COLA_check=zeros(length(input_signal),2);

frame_length=512;

for frame_shift=1:frame_length

    % Amplitude normalization for imposing unity COLA
    window=hanning (frame_length,'periodic');
    COLA_ratio=sum(window.*window)/frame_shift;
    window=window/sqrt(COLA_ratio);
    
    output_signal=zeros(length(input_signal),1);
    pin=0;     % current position in the input signal
    pout=0;   % current position in the output signal

    while pin+frame_length<length(input_signal)

        % Creating analysis frames
        analysis_frame = input_signal(pin+1:pin+frame_length).* window;
        
        % Leaving analysis frames untouched
        synthesis_frame=analysis_frame;
        
        % Weighted OverLap-Add (WOLA)
        output_signal(pout+1:pout+frame_length) = ...
            output_signal(pout+1:pout+frame_length)+synthesis_frame.*window;

        % Checking the COLA constaint for two values of the frame shift
        if (frame_shift==frame_length/2) 
            COLA_check(pout+1:pout+frame_length,1)=...
                COLA_check(pout+1:pout+frame_length,1)+window.*window;
        elseif (frame_shift==frame_length/4) 
            COLA_check(pout+1:pout+frame_length,2)=...
                COLA_check(pout+1:pout+frame_length,2)+window.*window;
        end;

        pin=pin+frame_shift;
        pout=pout+frame_shift;

    end;

    % Storing the output signal for frame shift = frame_length/2
        if (frame_shift==frame_length/2) 
            output_half=output_signal; 
        end;
    
    % Using the homemade |snr| function introduced in Chapter 2, 
    % and dropping the first and last frames, which are only partially
    % overlap-added.
    snr_values(frame_shift)=snr(input_signal(frame_length+1:end-frame_length),...
        output_signal(frame_length+1:end-frame_length),0);

end

plot((1:frame_length)/frame_length,snr_values);
xlabel('frame shift/frame length');
ylabel('SNR (dB)');

%%
% When the input signal is a unity constant, the output signal can indeed
% be close to unity... or not, depending the value of the frame shift.
% For the Hanning window, obviously, a frame shift of N/4 meets the COLA
% constraint, while N/2 leads to an important amplitude ripple.
% These values also depend on the weighting window (see Appendix 1 for a
% test on the square-root Hanning window)

plot(COLA_check(:,1)); hold on;
plot(COLA_check(:,2),'k--','linewidth',2); hold off;
xlabel('samples')
legend ('Frame shift = N/2','Frame shift = N/4')

%%
% A bad choice of the frame shift may lead to important perceptual
% degradation of the output audio signal. 

soundsc(output_half, Fs);
specgram(output_half,1024,Fs,256);


%% 1.3 STFT-based signal processing
% We now add an STFT/ISTFT step. The choice of the frame length N, of the
% type of weighting window and of the frame shift, depends on the type of
% processing applied to frequency bands. For frequency-selective
% processing, a large value of N is required, as the bandwidth of the DFT
% channels is inversely proportional to N. In this Section, we will
% implement a simple robotization effect on a speech signal, which is not a
% frequency-selective modification. We therefore set N to 256 samples (but
% many other values would match our goal). We choose a Hanning window, with
% a frame shift of 80 samples. This value does not meet the COLA
% constraint for the Hanning window, but the effect we will apply here does
% not attempt to maintain the integrity of the input signal.
% Let us process the 'speech.wav' file, which contains the sentence "Paint
% the circuits" sampled at 8 kHz. We first plot its spectrogram using the
% same frame length and frame shift as that of our STFT. The resulting
% plot is not especially pretty, but it shows exactly the data that we will
% process. The spectrogram reveals the harmonic structure of the signal,
% and shows that its pitch is a function of time.

frame_length=256;
frame_shift=80;
window=hanning (frame_length,'periodic');
COLA_ratio=sum(window.*window)/frame_shift;
window=window/sqrt(COLA_ratio);

[input_signal,Fs]=wavread('speech.wav');
specgram(input_signal,frame_length,Fs,window);
soundsc(input_signal,Fs);

%%
% When no processing of the intermediate DFTs is performed, the output
% signal is thus equal to the input signal.

output_signal=zeros(length(input_signal),1);
pin=0;pout=0;

while pin+frame_length<length(input_signal)

    % STFT
    analysis_frame = input_signal(pin+1:pin+frame_length).* window;
    dft_frame=fft(analysis_frame);
    
    % No processing
    
    % ISTFT
    synthesis_frame=ifft(dft_frame).*window;
	output_signal(pout+1:pout+frame_length) = ...
        output_signal(pout+1:pout+frame_length)+synthesis_frame;
    
    pin=pin+frame_shift;
    pout=pout+frame_shift;

end;

specgram(output_signal,frame_length,Fs,window);
soundsc(output_signal,Fs);

%% 
% It is now easy to modify the amplitudes of phases of the STFT,
% to produce audio effects. Let us test a simple robotization effect, for
% instance, by setting all phases to zero. This produces a disruptive
% perceptual effect, at the periodicity of the frame shift. When the shift
% is chosen small enough, and the effect is applied to speech it is
% perceived as artificial constant pitch.

output_signal=zeros(length(input_signal),1);
pin=0;pout=0;

while pin+frame_length<length(input_signal)

    % STFT
    analysis_frame = input_signal(pin+1:pin+frame_length).* window;
    dft=fft(analysis_frame);
    
    % Setting all phases to zero
    dft=abs(dft);
    
    % ISTFT
    synthesis_frame=ifft(dft).*window;
	output_signal(pout+1:pout+frame_length) = ...
        output_signal(pout+1:pout+frame_length)+synthesis_frame;
    
    pin=pin+frame_shift;
    pout=pout+frame_shift;

end;

% Zoom on a few ms of signal, before and after robotization
clf;
range=(400:1400);
subplot(211);
plot(range/Fs,input_signal(range))
subplot(212);
plot(range/Fs,output_signal(range))
xlabel('Time [ms]');

%%
%This simple effect appears on the spectrogram as a horizontal reshaping of
% the harmonics. It reveals the importance of the phases in the perception
% of a sound. 

clf
specgram(output_signal,frame_length,Fs,window);
soundsc(output_signal,Fs);

%% 2. Time scale modification
% In this Section we examine methods for modifying the duration of an input
% signal without changing its audio spectral characteristics (e.g., without
% changing the pitch of any of the instruments).
%
% The test sentence we use here contains several chunks, sampled at 44100
% Hz. It starts with 1000 Hz sine wave, followed by a speech exerpt. The
% sound of a (single) violin comes next, followed by an exerpt of a more
% complex polyphonic musical piece ('Time' by Pink Floyd).

[input_signal,Fs]=wavread('time_scaling.wav');

% (Again, for the spectrogram, we use the same frame length, frame shift, and
% weighting window as those we will use later in the phase vocoder)
frame_length=2048;
frame_shift=frame_length/4;
window=hanning (frame_length,'periodic');
COLA_ratio=sum(window.*window)/frame_shift;
window=window/sqrt(COLA_ratio);

soundsc(input_signal,Fs);
specgram(input_signal,frame_length,Fs,window);

%% 2.1 Interpolating the signal
% Increasing the length of the input signal by a factor of two, for
% instance, is easily obtained by interpolating the signal, using the
% |resample| function provided by MATLAB, and not telling the DAC about the
% sampling frequency change. This operation, however, also changes the
% frequency content of the signal: all frequencies are divided by two. 
% Similarly, speeding-up the signal by a factor two would multiply all
% frequencies by two. This is sometimes referred to as the "chipmunk
% effect".
%
% Notice incidentally the aliasing introduced by the imperfect low-pass
% filter used by |resample|. 

resampled_signal=resample(input_signal,2,1);

specgram(resampled_signal,frame_length,Fs,window);
soundsc(resampled_signal,Fs);

%% 2.2 Time-scale modification with Weighted OverLap-Add (WOLA)
% It is also possible to modify the time scale by decomposing it into
% ovelapping frames, changing the analysis frame shift into a different
% synthesis frame shift, and applying weighted overlap-add (WOLA) to the
% resulting synthesis frames. While this technique works well for
% unstructured signals, its application to harmonic signals produces
% unpleasant audio artifacts at the frequency of the synthesis frame rate,
% due to the loss of synchronicity between excerpts of the same partials in
% overlapping synthesis frames. 
% This somehow reminds us of the robotization effect obtained in Section
% 1.3. 

frame_length=2048;
synthesis_frame_shift=frame_length/4;
window=hanning (frame_length,'periodic');
COLA_ratio=sum(window.*window)/synthesis_frame_shift;
window=window/sqrt(COLA_ratio);

time_scaling_ratio=2.85;
analysis_frame_shift=round(synthesis_frame_shift/time_scaling_ratio);

pin=0;pout=0;
output_signal=zeros(time_scaling_ratio*length(input_signal),1);

while (pin+frame_length<length(input_signal)) ...
        && (pout+frame_length<length(output_signal))

    analysis_frame = input_signal(pin+1:pin+frame_length).* window;
    synthesis_frame = analysis_frame.* window;
	output_signal(pout+1:pout+frame_length) = ...
        output_signal(pout+1:pout+frame_length)+synthesis_frame;

    pin=pin+analysis_frame_shift;
    pout=pout+synthesis_frame_shift;
    
    % Saving frames for later use
    if (pin==2*analysis_frame_shift) % 3rd frame
        frame_3=synthesis_frame;
    elseif (pin==3*analysis_frame_shift) % 4th frame
        frame_4=synthesis_frame;
    end;
    
end;

% Plot two overlapping synthesis frames and show the resulting output signal
clf;
ax(1)=subplot(211);
range=(2*synthesis_frame_shift:2*synthesis_frame_shift+frame_length-1);
plot(range/Fs,frame_3)
hold on;
range=(3*synthesis_frame_shift:3*synthesis_frame_shift+frame_length-1);
plot(range/Fs,frame_4,'r')
hold off
ax(2)=subplot(212);
range=(2*synthesis_frame_shift:3*synthesis_frame_shift+frame_length-1);
plot(range/Fs,output_signal(range))
xlabel('Time [s]');
linkaxes(ax,'x');
set(gca,'xlim',[0.045 0.06]);

soundsc(output_signal,Fs);

%% 2.3 Time-scale modification with the phase vocoder
% Let us now modify the time scale of the input signal without
% affecting its frequency content, by using a phase vocoder, which is a
% STFT-based signal processing system with specific hypotheses on the STFT.
% Time-scaling is again achieved by modifying the analysis frame shift without
% changing the synthesis frame shift, but the STFT has to be modified, so as to
% avoid creating phase mismatches between overlapped synthesis frames.
%
% We use the periodic Hanning window, which provides a low spectral leakage 
% while its main spectral lobe width is limited to 4*Fs/N. The choice of
% the frame length N is dictated by a maximum bandwidth constraint for the
% DFT sub-band filters, since the phase vocoder is based on the hypothesis
% that each DFT bin is mainly influenced by a single partial, i.e. a
% sinusoïdal component with constant frequency. For a Hanning window, each 
% sinusoidal component will drive 4 DFT bins, the width of each is Fs/N.
% Clearly, the higher the value of N, the better the selectivity of the
% sub-band filters. On the other hand, a high value of N tends to break the
% stationarity hypothesis on the analysis frame (i.e., the frequency of
% partials will no longer be constant). We set N to 2048 samples here,
% which imposes the width of each DFT bin to 21 Hz, and the bandwidth of
% the DFT sub-band filters to 86 Hz. We set the frame shift to N/4=512
% samples, which meets the COLA constraint for the Hanning window. What is
% more, this choice obeys the Shannon theorem for the output of each DTF
% bin, seen as the output of a sub-band filter.

NFFT=2048;
frame_length=NFFT;
synthesis_frame_shift=frame_length/4;
window=hanning (frame_length,'periodic');
COLA_ratio=sum(window.*window)/synthesis_frame_shift;
window=window/sqrt(COLA_ratio);

time_scaling_ratio=2.85;
analysis_frame_shift=round(synthesis_frame_shift/time_scaling_ratio);
% Central frequency of each DFT channel
DFT_bin_freqs=((0:NFFT/2-1)*2*pi/NFFT)';

pin=0;pout=0;

output_signal=zeros(time_scaling_ratio*length(input_signal),1);
last_analysis_phase=zeros(NFFT/2,1);
last_synthesis_phase=zeros(NFFT/2,1);

while (pin+frame_length<length(input_signal)) ...
        && (pout+frame_length<length(output_signal))

    % STFT
    analysis_frame = input_signal(pin+1:pin+frame_length).* window;
    dft=fft(fftshift(analysis_frame));
    dft=dft(1:NFFT/2);

    % PHASE MODIFICATION
    % Find phase for each bin, compute by how much it increased since last
    % frame. 
    this_analysis_phase = angle(dft); 
    delta_phase = this_analysis_phase - last_analysis_phase;
    phase_increment=delta_phase-analysis_frame_shift*DFT_bin_freqs;
    
    % Estimate the frequency of the main partial for each bin
    principal_determination=mod(phase_increment+pi,2*pi)-pi;
    partials_freq=principal_determination/analysis_frame_shift+DFT_bin_freqs;
    
    % Update the phase in each bin
    this_synthesis_phase=last_synthesis_phase+synthesis_frame_shift*partials_freq;
    
    % Compute DFT of the synthesis frame
    dft= abs(dft).* exp(j*this_synthesis_phase);
    
    % Remember phases
    last_analysis_phase=this_analysis_phase;
    last_synthesis_phase=this_synthesis_phase;
    
    % ISTFT
    dft(NFFT/2+2:NFFT)=fliplr(dft(2:NFFT/2)');
    synthesis_frame = fftshift(real(ifft(dft))).* window;
	output_signal(pout+1:pout+frame_length) = ...
        output_signal(pout+1:pout+frame_length)+synthesis_frame;

    pin=pin+analysis_frame_shift;
    pout=pout+synthesis_frame_shift;

    % Saving the estimated frequency of partials for later use
    if (pin==2*analysis_frame_shift) % 3rd frame
        partials_freq_3=partials_freq;
    end;
    % Saving frames for later use
    if (pin==2*analysis_frame_shift) % 3rd frame
        frame_3=synthesis_frame;
    elseif (pin==3*analysis_frame_shift) % 4th frame
        frame_4=synthesis_frame;
    end;
    
end;

clf;
specgram(output_signal(1:end-5000),frame_length,Fs,window);

%% 
% The phasing effect we had observed in the previous method has now disappeared.

% Plot two overlapping synthesis frames and show the resulting output signal
clf;
ax(1)=subplot(211);
range=(2*synthesis_frame_shift:2*synthesis_frame_shift+frame_length-1);
plot(range/Fs,frame_3)
hold on;
range=(3*synthesis_frame_shift:3*synthesis_frame_shift+frame_length-1);
plot(range/Fs,frame_4,'r')
hold off
ax(2)=subplot(212);
range=(2*synthesis_frame_shift:3*synthesis_frame_shift+frame_length-1);
plot(range/Fs,output_signal(range))
xlabel('Time [s]');
linkaxes(ax,'x');
set(gca,'xlim',[0.045 0.06]);

soundsc(output_signal,Fs);

%%
% Let us check how the reduced frequency of the initial sinusoid (i.e.
% 1000*2*pi/Fs = 0.142 rad/s) has been estimated. 
% The first spectral line of the sinusoid mainly influences DTF bins around
% index 1000/(44100/2048)=46.43. As observed on the plot, bins 42
% to 51 (which correspond to the main spectral lobe) have indeed measured
% the correct frequency. 

subplot(2,1,1);
plot((0:NFFT/2-1),partials_freq_3);
set(gca,'xlim',[25 67]);
ylabel('Frequency (rad/s)');

subplot(2,1,2);
third_frame=input_signal(2*analysis_frame_shift:2*analysis_frame_shift+frame_length-1);
dft=fft(third_frame.* window);
dft=20*log10(abs(dft(1:NFFT/2)));
plot((0:NFFT/2-1),dft);
ylabel('Amplitude (dB)');xlabel('DFT bin');
set(gca,'xlim',[25 67]);

%%
% However, for more complex sounds (as in the last frame of our test
% signal), the estimated partial frequency in neighboring bins vary a lot. 

subplot(2,1,1);
plot((0:NFFT/2-1),partials_freq);
set(gca,'xlim',[25 67]);
ylabel('Frequency (rad/s)');

subplot(2,1,2);
dft=fft(analysis_frame);
dft=20*log10(abs(dft(1:NFFT/2)));
plot((0:NFFT/2-1),dft);
ylabel('Amplitude (dB)');xlabel('DFT bin');
set(gca,'xlim',[25 67]);

%%
% As a result, although the spectrogram of the time-scaled signal looks similar to
% that of the original signal (except of course for the time axis), it
% exhibits significant phasiness and transient smearing.

soundsc(output_signal,Fs);

%% 2.4 The phase-locked vocoder
% We have seen in the previous Section that, in the case of a sinusoidal
% input signal, the estimated frequency of the main partial in several 
% neighboring DFT bins (around the frequency of the sinousoid) is somehow
% constant. One of the reasons of the phasiness in the the previous verison
% of the phase vocoder comes from the fact that it does not enforce this
% effect. If other small amplitude partials are added the sinusoid, each
% DFT bin computes its own partial frequency, so that these estimations in
% neighboring bins will only coincide by chance. The phase-locked vocoder
% changes this, by at locking the estimation of partial frequencies in bins
% surrounding spectral peaks. This is termed as vertical phase-locking.
% In the following lines, we give an implementation of the identity
% phase-locking scheme proposed in (Laroche and Dolson 1999).
%
% *MATLAB function involved:*
% 
% * |function [peaks,regions] = findpeaks(x)| returns the indexes of the
% local maxima (|peaks|) and of the minima (|regions|) in array |x|.

pin=0;pout=0;

output_signal=zeros(time_scaling_ratio*length(input_signal),1);
last_analysis_phase=zeros(NFFT/2,1);
last_synthesis_phase=zeros(NFFT/2,1);

while (pin+frame_length<length(input_signal)) ...
    && (pout+frame_length<length(output_signal))

    % STFT
    analysis_frame = input_signal(pin+1:pin+frame_length).* window;
    dft=fft(fftshift(analysis_frame));
    dft=dft(1:NFFT/2);
    
    % PHASE MODIFICATION
    % Find phase for each bin, compute by how much it increased since last
    % frame. 
    this_analysis_phase = angle(dft); 
    delta_phase = this_analysis_phase - last_analysis_phase;
    phase_increment=delta_phase-analysis_frame_shift*DFT_bin_freqs;

    % Find peaks in DFT.
    peaks = findpeaks(abs(dft));

    % Estimate the frequency of the main partial for each bin
    principal_determination(peaks)=mod(phase_increment(peaks)+pi,2*pi)-pi;
    partials_freq(peaks)=principal_determination(peaks)/...
        analysis_frame_shift+DFT_bin_freqs(peaks);

    % Find regions of influence around peaks
    regions = round(0.5*(peaks(1:end-1)+peaks(2:end)));
    regions = [1;regions;NFFT/2];

    % Set the frequency of partials in regions of influence to that of
    % the peak (this is not strictly needed; it is used for subsequent plots)
    for i=1:length(peaks)
        partials_freq(regions(i):regions(i+1)) = partials_freq(peaks(i));
    end

    % Update the phase for each peak
    this_synthesis_phase(peaks)=last_synthesis_phase(peaks)+...
        synthesis_frame_shift*partials_freq(peaks);

    % force identity phase locking in regions of influence
    for i=1:length(peaks)
        this_synthesis_phase(regions(i):regions(i+1)) = this_synthesis_phase(peaks(i))...
            + this_analysis_phase(regions(i):regions(i+1))-this_analysis_phase(peaks(i));
    end

    % Compute DFT of the synthesis frame
    dft= abs(dft).* exp(j*this_synthesis_phase);
    % Remember phases
    last_analysis_phase=this_analysis_phase;
    last_synthesis_phase=this_synthesis_phase;

    % ISTFT
    dft(NFFT/2+2:NFFT)=fliplr(dft(2:NFFT/2)');
    synthesis_frame = fftshift(real(ifft(dft))).* window;
	output_signal(pout+1:pout+frame_length) = ...
        output_signal(pout+1:pout+frame_length)+synthesis_frame;

    pin=pin+analysis_frame_shift;
    pout=pout+synthesis_frame_shift;

    % Saving the estimated frequency of partials for later use
    if (pin==2*analysis_frame_shift) % third frame
        partials_freq_3=partials_freq;
    end
    
end;

clf;
specgram(output_signal(1:end-5000),frame_length,Fs,window);

%%
% As a result of vertical phase-locking, a much larger number of bins
% around the main spectral line of our sinusoid have been assigned the same
% partial frequency.

subplot(2,1,1);
plot((0:NFFT/2-1),partials_freq_3);
ylabel('Frequency (rad/s)');
set(gca,'xlim',[25 67]);

subplot(2,1,2);
third_frame=input_signal(2*analysis_frame_shift:2*analysis_frame_shift+frame_length-1);
dft=fft(third_frame.* window);
dft=20*log10(abs(dft(1:NFFT/2)));
plot((0:NFFT/2-1),dft);
ylabel('Amplitude (dB)');xlabel('DFT bin');
set(gca,'xlim',[25 67]);

%%
% More importantly, this property has been enforced in more complex sounds,
% as shown in the last frame of our test signal.

subplot(2,1,1);
plot((0:NFFT/2-1),partials_freq);
ylabel('Frequency (rad/s)');
set(gca,'xlim',[25 67]);

subplot(2,1,2);
dft=fft(analysis_frame);
dft=20*log10(abs(dft(1:NFFT/2)));
plot((0:NFFT/2-1),dft);
ylabel('Amplitude (dB)');xlabel('DFT bin');
set(gca,'xlim',[25 67]);

%%
% The result is improves, although it is still far from natural. 

soundsc(output_signal,Fs);

%%
% It is possible to attenuate transient smearing by using separate analysis
% windows for stationary and transient sounds (as done in the MPEG audio
% coders). We do not examine this option here.

%% 3. Pitch modification
%
% The simplest way to produce a modification of the frequency axis of a
% signal is simply to resample it, as already seen in Section 2.1. This
% method, however, also produces time modification.
% In this Section we focus on methods for modifying the pitch of an input
% signal without changing its duration.  

%% 3.1 Time-scale modification and resampling 
%
% In order to avoid the time-scaling effect of a simple resampling-based
% approach, a compensatory time-scale modification can be applied, using a
% phase vocoder.  
% Let us multiply the pitch by a factor 0.7, for instance. We first
% multiply the duration of the input signal by 0.7, and then resample it by
% 10/7, so as to recover the original number of samples.
%
% *MATLAB function involved:*
% 
% * |[output_signal] = phase_locked_vocoder(input_signal,time_scaling_ratio)|
% is a simple implementation of the phase locked vocoder described in:
% Laroche J, Dolson M (1999) Improved Phase Vocoder Time-Scale Modification
% of Audio. IEEE Trans. Speech and Audio Processing, 3, pp 323–332

[input_signal,Fs]=wavread('time_scaling.wav');

time_scaling_ratio=0.7;
resampling_ratio=2;

output_signal=phase_locked_vocoder(input_signal,time_scaling_ratio);
pitch_modified_signal=resample(output_signal,10,7);

soundsc(input_signal,Fs);
soundsc(pitch_modified_signal,Fs);

%% 
% The time scale of the output signal is now identical to that of the input
% signal. Notice again the aliasing introduced by the imperfect low-pass
% filter used by |resample|. 

clf
specgram(pitch_modified_signal(1:end-5000),frame_length,Fs,window);

%% 3.2 The STFT-based approach
%
% It is also possible, and much easier, to perform pitch modification in
% the frequency domain by translating partials to new frequencies, i.e.
% inside the phase vocoder. This technique is much more flexible, as it
% allows for non-linear modification of the frequency axis.

NFFT=2048;
frame_length=NFFT;
frame_shift=frame_length/4;
window=hanning (frame_length,'periodic');
COLA_ratio=sum(window.*window)/frame_shift;
window=window/sqrt(COLA_ratio);

pitch_scaling_ratio=0.7;

% Central frequency of each DFT channel
DFT_bin_freqs=((0:NFFT/2-1)*2*pi/NFFT)';

pin=0;pout=0;

output_signal=zeros(length(input_signal),1);
accumulated_rotation_angles=zeros(1,NFFT/2);

while pin+frame_length<length(input_signal)

    % STFT
    analysis_frame = input_signal(pin+1:pin+frame_length).* window;
    dft=fft(fftshift(analysis_frame));
    dft=dft(1:NFFT/2);
    
    % Find peaks in DFT.
    peaks = findpeaks(abs(dft));

    % Find regions of influence around peaks
    regions = round(0.5*(peaks(1:end-1)+peaks(2:end)));
    regions = [1;regions;NFFT/2];

    % Move each peak in frequency according to modification factor.
    %
    % NB: This is the simplest possible version: the peak frequencies are
    % not interpolated for better precision, and the shifting is rounded to
    % an integer number of bins. This will produce modulation artifacts on
    % sweeping sinusoids, and phasing artifacts on audio signals such as
    % speech. Moreover, when shifting the pitch down, peaks are truncated
    % at DC, whereas they should be reflected with a conjugation sign.

    modified_dft = zeros(size(dft));
    for u=1:length(peaks)
        % Locate old and new bin.
        old_bin = peaks(u)-1;
        new_bin = round(pitch_scaling_ratio*old_bin);
        % Be sure to  stay within 0-NFFT/2 when shifting/copying the peak
        % bins
        if(new_bin-old_bin+regions(u) >= NFFT/2) break; end;
        if(new_bin-old_bin+regions(u+1) >= NFFT/2) 
            regions(u+1) = NFFT/2 - new_bin + old_bin; 
        end;
        if(new_bin-old_bin+regions(u) <= 0) 
            regions(u) = 1 - new_bin + old_bin; 
        end;
        % Compute the rotation angle required, which has to be cumulated from
        % frame to frame
        rotation_angles = accumulated_rotation_angles(old_bin+1) ...
            + 2*pi*frame_shift*(new_bin-old_bin)/NFFT;
        % Overlap/add the bins around the peak, changing the phases
        % accordingly 
        modified_dft(new_bin-old_bin+(regions(u):regions(u+1))) = ...
            modified_dft(new_bin-old_bin+(regions(u):regions(u+1))) + ...
            dft(regions(u):regions(u+1)) * exp(j*rotation_angles);
        accumulated_rotation_angles((regions(u):regions(u+1))) = ...
            rotation_angles;
    end


    % ISTFT
    modified_dft(NFFT/2+2:NFFT)=fliplr(modified_dft(2:NFFT/2)');
    synthesis_frame = fftshift(real(ifft(modified_dft))).* window;
	output_signal(pout+1:pout+frame_length) = ...
        output_signal(pout+1:pout+frame_length)+synthesis_frame;

    pin=pin+frame_shift;
    pout=pout+frame_shift;

end;

soundsc(output_signal,Fs);

%% 
% The time scale of the output signal is still identical to that of the
% input signal, and the frequency content above |Fs*time_scaling_ratio| is
% set to zero.

clf
specgram(output_signal(1:end-5000),frame_length,Fs,window);

%% Appendix 1 - The square-root Hanning window
% Although we have systematically used the Hanning window in the above
% proof-of-concept, it is worth mentioning that the square-root of the
% ususal windows can also be used in STFT-based signal processing.
% The square-root Hanning window is a good example. Let us compare its
% waveform and spectral features to those of the usual wieghting windows.
% It appears that the spectral leakage (i.e. the level of sidelobes) is
% higher than with the Hanning window, which explains why the latter is
% preferred in phase vocoders. 

N=100;
wvtool(boxcar(N),sqrt(hanning(N)),hanning(N),blackman(N));

%% COLA constraint
% It is easy to check that the sqrt(Hanning) window does meet the COLA
% constraint for N/2. Besides, this window, with this specific frame shift,
% is used in the Modulated Lapped Transform (MLT) in perceptpual audio
% coding (see Chapter 3).

Fs=8000;
input_signal=chirp((0:Fs)/Fs,0,1,2000)';
specgram(input_signal,1024,Fs,256);

COLA_check=zeros(length(input_signal),2);

frame_length=512;

for frame_shift=1:frame_length

    % Amplitude normalization for imposing unity COLA
    window=sqrt(hanning (frame_length,'periodic'));
    COLA_ratio=sum(window.*window)/frame_shift;
    window=window/sqrt(COLA_ratio);
    
    output_signal=zeros(length(input_signal),1);
    pin=0;     % current position in the input signal
    pout=0;   % current position in the output signal

    while pin+frame_length<length(input_signal)

        % Creating analysis frames
        analysis_frame = input_signal(pin+1:pin+frame_length).* window;
        
        % Leaving analysis frames untouched
        synthesis_frame=analysis_frame;
        
        % Weighted OverLap-Add (WOLA)
        output_signal(pout+1:pout+frame_length) = ...
            output_signal(pout+1:pout+frame_length)+synthesis_frame.*window;

        % Checking the COLA constaint for two values of the frame shift
        if (frame_shift==frame_length/2) 
            COLA_check(pout+1:pout+frame_length,1)=...
                COLA_check(pout+1:pout+frame_length,1)+window.*window;
        elseif (frame_shift==frame_length/4) 
            COLA_check(pout+1:pout+frame_length,2)=...
                COLA_check(pout+1:pout+frame_length,2)+window.*window;
        end;

        pin=pin+frame_shift;
        pout=pout+frame_shift;

    end;
    
    % Using the homemade |snr| function introduced in Chapter 2, 
    % and dropping the first and last frames, which are only partially
    % overlap-added.
    snr_values(frame_shift)=snr(input_signal(frame_length+1:end-frame_length),...
        output_signal(frame_length+1:end-frame_length),0);

end

plot((1:frame_length)/frame_length,snr_values);
xlabel('frame shift/frame length');
ylabel('SNR (dB)');

