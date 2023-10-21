%% 3. 32-Channel Pseudo-QMF filter bank
hn=PQMF32_prototype;

%%

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

[input_signal,Fs]=audioread('violin.wav');
N = length(input_signal)
slen = N/Fs

[y,Fs]=audioread('01 David Bowie - Blackstar.wav');
cut = y(1:slen*Fs,:);
mono = sum(cut,2)/2;

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



clf;
figure
specgram(G3_output,1024,Fs,256);
soundsc(G3_output,Fs);

figure
specgram(output_signal,1024,Fs,256);
soundsc(output_signal,Fs);


error=output_signal(512:end)-mono(1:end-511);

[signal_psd,w]=periodogram(mono(11001:12024),hamming(1024));
[error_psd,w]=periodogram(error(11001:12024),hamming(1024));
figure
plot(w/pi*22050,10*log10(signal_psd)); hold on;
plot(w/pi*22050,10*log10(error_psd),'r','linewidth',2); hold off;
set(gca,'xlim',[0 22050]);
legend('Signal PSD', 'Error PSD');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');

snr_PQMF=snr(mono(1:end-511),output_signal(512:end),0)

