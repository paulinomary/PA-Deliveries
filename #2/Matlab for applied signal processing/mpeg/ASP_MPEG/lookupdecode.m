function [indexmessage, codemessage] = lookupdecode(codemessage, decoder, n, option)

%lookupdecode   Lookup table decoding (for internal use only).
%
%   [indexmessage, codemessage] = 
%       lookupdecode(codemessage, decoder[, n[, option]]) decodes a message coded
%   with a prefix code into an index message using the lookup decoder 'decoder'.
%   Decodes at most 'n' symbols.  If 'n' is not specified, it is taken to be Inf,
%   i.e., all possible symbols are decoded.
%
%   Hashmarks are shown during decoding if the option 'verbose' is passed.
%
%   See also decode, treedecode, lookupdecoder, treedecoder.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if ~isscalar(decoder) | ~isstruct(decoder)
    error('decoder must be a scalar struct');
elseif isempty(strmatch('type', fieldnames(decoder)))
    error('decoder must have a type field');
elseif ~strcmp(decoder.type, 'lookup')
    error('only a lookup decoder may be passed to this function');
elseif isempty(strmatch('n_symbols', fieldnames(decoder)))
    error('decoder must have a n_symbols field');
elseif decoder.n_symbols < 2
    error('decoder must correspond to a code with two possible codewords');
end

if nargin < 3
    n = Inf;
end

if nargin < 4
    verbose = false;
else
    if ~isa(option, 'char')
        error('option must be a string');
    end
    verbose = strcmp(option, 'verbose');
end

% Decode:

if(n ~= Inf)
    indexmessage = zeros(1, n);
else
    % Initialize with some default size:
    indexmessage = zeros(1, 1000);
end

nbits = decoder.nbits;
table = decoder.table;
slowdecoder = decoder.slowdecoder;

if n == Inf
    hash = 50;
else
    hash = ceil(n / 79);
end

s = 1;
b = 1;
while(s <= n & b + nbits - 1 <= length(codemessage))
    prevs = s;
    cm = bin2dec(codemessage(b : b + nbits - 1));
    % mms 2002/5/22: Check whether cm is decodeable.  If it is not, abort.
    if ~table(cm + 1).isvalid
        error('code message is not decodeable (code is redundant).');
    end
    if(length(table(cm + 1).indices) == 0)
        % Oppps...  Longish codeword...  Try to decode it with normal decode.
        [im, rem] = feval(slowdecoder.decode, ...
            codemessage(b : end), slowdecoder, 1);
        if(length(im) == 0)
            break;
        else
            m = length(im);
            if(length(indexmessage) < s + m - 1)
                indexmessage(max([s + m - 1  2 * length(indexmessage)])) = 0;
            end
            indexmessage(s : s + m - 1) = im;
            s = s + m;
            b = b + length(decoder.code{im});
        end
    else
        m = length(table(cm + 1).indices);
        if(length(indexmessage) < s + m - 1)
            indexmessage(max([s + m - 1  2 * length(indexmessage)])) = 0;
        end
        indexmessage(s : s + m - 1) = table(cm + 1).indices;
        s = s + m;
        b = b + nbits - table(cm + 1).remaining;
    end
    
    if verbose
        if prevs < s
            for i = prevs + 1 : s
                if mod(i, hash) == 0
                    fprintf('#');
                    if mod(floor(i / hash), 79) == 0
                        fprintf('\n');
                    end
                end
            end
        end
    end
end

if verbose
    fprintf('\n');
end

codemessage = codemessage(b : end);

% Try to decode remaining symbols, if there are symbols missing:

if(n >= s) 
    [im, codemessage] = feval(slowdecoder.decode, ...
        codemessage, slowdecoder, n - s + 1);

    m = length(im);
    if(length(indexmessage) < s + m - 1)
        indexmessage(max([s + m - 1  2 * length(indexmessage)])) = 0;
    end
    indexmessage(s : s + m - 1) = im;
    s = s + m;
end
        
if(s > n + 1)
    cm = encode(indexmessage(n + 1 : s - 1), decoder.code);
    codemessage = [cm codemessage];
    indexmessage = indexmessage(1 : n);
else
    indexmessage = indexmessage(1 : s - 1);
end