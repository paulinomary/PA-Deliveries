function vdct = iquantizedct8diff(levels, step)

%iquantizedct8diff   Values from uniform quantization levels for differencial DCT
%                    coeficients.
%
%   vdct = iquantizedct8diff(levels, step) returns the recovered values from 
%   quantization levels 'levels' obtained using uniform quantization with
%   quantization steps 'step'.  'levels' is assumed to be a 8x8 uint8 matrix.
%   'step' may either be a scalar or an array with the same size as 'levels'.  
%   The quantization levels are assumed to originate from the coeficients of a 
%   two-dimensional DCT calculated over 8x8 blocks of diferences between uint8 
%   images.  The quantization levels are also assumed to be offset by 128.
%
%   See quantizedct8diff, quantizedct8, and iquantizedct8.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

levels = double(levels) - 128;

vdct = iquantize(levels, step);
