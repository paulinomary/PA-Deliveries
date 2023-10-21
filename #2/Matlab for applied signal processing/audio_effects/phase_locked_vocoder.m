function [output_signal] = phase_locked_vocoder(input_signal,time_scaling_ratio)

% |[output_signal] = phase_locked_vocoder(input_signal,time_scaling_ratio)|
% is a simple implementation of the phase locked vocoder described in:
% Laroche J, Dolson M (1999) Improved Phase Vocoder Time-Scale Modification
% of Audio. IEEE Trans. Speech and Audio Processing, 3, pp 323–332

% T DUTOIT, 09/04/2008

NFFT=2048;
frame_length=NFFT;
synthesis_frame_shift=frame_length/4;
window=hanning (frame_length,'periodic');
COLA_ratio=sum(window.*window)/synthesis_frame_shift;
window=window/sqrt(COLA_ratio);

analysis_frame_shift=round(synthesis_frame_shift/time_scaling_ratio);
% Central frequency of each DFT channel
DFT_bin_freqs=((0:NFFT/2-1)*2*pi/NFFT)';

pin=0;pout=0;

output_signal=zeros(time_scaling_ratio*length(input_signal),1);
this_analysis_phase=zeros(NFFT/2,1);
principal_determination=zeros(NFFT/2,1);
partials_freq=zeros(NFFT/2,1);
this_synthesis_phase=zeros(NFFT/2,1);
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

end;
