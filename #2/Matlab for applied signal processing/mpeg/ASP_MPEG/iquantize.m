function v = iquantize(levels, step)

%iquantize   Values from uniform quantization levels.
%
%   v = iquantize(levels, step) returns the recovered values from quantization
%   levels 'levels' obtained using uniform quantization with quantization steps 
%   'step'.  'levels' may have any dimensions and size.  'step' may either be a
%   scalar or an array with the same size as 'levels'.
%
%   See quantize.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

v = levels .* step;