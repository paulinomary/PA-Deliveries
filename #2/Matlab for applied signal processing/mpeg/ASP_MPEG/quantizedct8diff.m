function levels = quantizedct8diff(vdct, step)

%quantizedct8diff   Uniform quantization levels from differencial DCT coeficients.
%
%   levels = quantizedct8diff(vdct, step) returns quantization levels by 
%   uniformly quantizing values in 'vdct' according to steps in 'step'.  
%   Quantization levels are saturated to the range [0 255], but only after the
%   quantized two-dimensional DCT coeficients have been offset by 128.
%
%   'vdct' must be a 8x8 matrix resulting from the diference between to uint8 
%   images.  'step' may either be a scalar or an array with the same size as 
%   'vdct'.
%
%   See iquantizedct8diff, quantizedct8, and iquantizedct8.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

levels = quantize(vdct, step, -128, 127);

levels = uint8(levels + 128);
