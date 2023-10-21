function [indexmessage, codemessage] = treedecode(codemessage, decoder, n, option)

%treedecode   Tree decoding (for internal use only).
%
%   [indexmessage, codemessage] = treedecode(codemessage, decoder[, n[, option]])
%   decodes a message coded with a prefix code into an index message using the 
%   tree decoder 'decoder'.  Decodes at most 'n' symbols.  If 'n' is not
%   specified, it is taken to be Inf, i.e., all possible symbols are decoded.
%
%   Hashmarks are shown during decoding if the option 'verbose' is passed.
%
%   See also decode, lookupdecode, lookupdecoder, treedecoder.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if ~isscalar(decoder) | ~isstruct(decoder)
    error('decoder must be a scalar struct');
elseif isempty(strmatch('type', fieldnames(decoder)))
    error('decoder must have a type field');
elseif ~strcmp(decoder.type, 'tree')
    error('only a tree decoder may be passed to this function');
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

tree = decoder.tree;

if(n ~= Inf)
    indexmessage = zeros(1, n);
else
    % Initialize with some default size:
    indexmessage = zeros(1, 1000);
end

if n == Inf
    hash = 5;
else
    hash = ceil(n / 79);
end

node = 1;
s = 0;
b = 1;
lb = 1;
while s ~= n & b <= length(codemessage)
    % Advance:
    node = tree(node).branches(logical(codemessage(b)) + 1);
    b = b + 1;

    if(node == 0)
        error('code message is not decodeable (code is redundant).');
    end
    % Check if at a leave of the tree:
    if(tree(node).index ~= 0)
        if verbose
            if mod(s, hash) == 0
                fprintf('#');
                if mod(floor(s / hash + 1), 79) == 0
                    fprintf('\n');
                end
            end
        end
        % Since it is a terminal node, store index:
        s = s + 1;
        indexmessage(s) = tree(node).index;
        node = 1;
        lb = b;
    end
end

if verbose
    fprintf('\n');
end

codemessage = codemessage(lb : end);
indexmessage = indexmessage(1 : s);