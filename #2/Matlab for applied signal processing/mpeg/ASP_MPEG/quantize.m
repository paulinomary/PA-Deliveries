function levels = quantize(v, step, min, max)

%quantize   Uniform quantization levels from values.
%
%   levels = quantize(v, step[, min[, max]]) returns quantization levels by 
%   uniformly quantizing values in 'v' according to steps in 'step'.  
%   Quantization levels are saturated to the range ['min' 'max'].  'min' has 
%   default value -Inf and 'max' has default value 'Inf'.
%
%   'v' may have any dimensions and size.  'step', 'min' and 'max' may either
%   be scalars or arrays with the same size as 'v'.
%
%   See iquantize.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if nargin < 2 | 4 < nargin
    error('wrong number of arguments');
elseif nargin < 4
    if nargin < 3
        min = -Inf;
    end
    max = Inf;
end

levels = limit(round(v ./ step), min, max);