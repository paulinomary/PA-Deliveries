function vdct = iquantizedct8(levels, step)

%iquantizedct8   Values from uniform quantization levels for DCT coeficients.
%
%   vdct = iquantizedct8(levels, step) returns the recovered values from 
%   quantization levels 'levels' obtained using uniform quantization with
%   quantization steps 'step'.  'levels' is assumed to be a 8x8 uint8 matrix. 
%   'step' may either be a scalar or an array with the same size as 'levels'.  
%   The quantization levels are assumed to originate from the coeficients of a 
%   two-dimensional DCT calculated over uint8 8x8 blocks.  The quantization levels are
%   also assumed to be offset by 128 in the case of the AC coeficients.
%
%   See quantizedct8, quantizedct8diff, and iquantizedct8diff.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

offset = ones(8) * 128;
offset(1, 1) = 0;

levels = double(levels) - offset;

vdct = iquantize(levels, step);

