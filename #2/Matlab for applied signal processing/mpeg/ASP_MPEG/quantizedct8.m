function levels = quantizedct8(vdct, step)

%quantizedct8   Uniform quantization levels from DCT coeficients.
%
%   levels = quantizedct8(vdct, step) returns quantization levels by 
%   uniformly quantizing values in 'vdct' according to steps in 'step'.  
%   Quantization levels are saturated to the range [0 255], but only after the
%   quantized AC two-dimensional DCT coeficients have been offset by 128.
%
%   'vdct' must be a double 8x8 matrix.  'step' may either be a scalar or
%   an array with the same size as 'vdct'.
%
%   See iquantizedct8, quantizedct8diff, and iquantizedct8diff.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

min = ones(8) * -128;
min(1, 1) = 0;

max = ones(8) * 127;
max(1, 1) = 255;

levels = quantize(vdct, step, min, max);

offset = ones(8) * 128;
offset(1, 1) = 0;

levels = uint8(levels + offset);
