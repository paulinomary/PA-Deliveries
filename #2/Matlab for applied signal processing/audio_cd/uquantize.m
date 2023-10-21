function signal_quantized = uquantize(signal,N,saturation)

% function signal_quantized = uquantize(signal,N,saturation) quantizes
% signal uniformly on N bits between -saturation and +saturation, using a
% mid-rise quantizer and triangular dither with a range of twice the
% quantization step. Clipping is performed when saturation is reached.
% 
% T DUTOIT, 13:49 12/03/2007

q=2*saturation/2^N;
dither=(rand(size(signal))-0.5)*q+(rand(size(signal))-0.5)*q;
signal_dithered=signal+dither;
signal_quantized=floor(signal_dithered/q)*q+q/2;
% clipping the quantized signal to [saturation-q/2,-saturation+q/2]
clipping=saturation-q/2;
signal_quantized=min(signal_quantized,clipping);
signal_quantized=max(signal_quantized,-clipping);

