function decoder = lookupdecoder(code, nbits, option)

%lookupdecoder - Builds a lookup decoder for a prefix code.
%
%   decoder = lookupdecoder(code[, nbits[, option]]) returns the lookup decoder
%   of code 'code'.  Assumes the code is prefix.  Useful for passing to the 
%   function 'decode' so that the lookup decoder does not has to be generated at 
%   each decoding.  'nbits' is by default 8, and corresponds to the number of 
%   bits to use during lookup.  The larger it is, the faster the decoder
%   (however, there is a memory tradeof, since the required memory grows 
%   exponencialy with 'nbits').
%
%   Code 'code' must be a cell of char (row) arrays.  Each element of this array
%   is the binary codeword.
%
%   Hashmarks are shown during the construction of the lookup table if the option 
%   'verbose' is passed.
%
%   Notice that a lookup decoder reverts to a tree decoder whenever the number of
%   available bits is insuficient for lookup.  A tree decoder is also used for
%   building the lookup table entries.  Hence, building a lookup decoder is a
%   slow operation, which preferably should be done only once.
%
%   See also decode, lookupdecode, treedecode, treedecoder.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if ~iscell(code) | ~isrow(code)
    error('code must be a row cell vector');
end

if(nargin < 2)
    nbits = 8;
elseif ~isscalar(nbits) | ~isinteger(nbits) | nbits < 0
    error('nbits must be a scalar positive integer')
end

if nargin < 3
    verbose = false;
else
    if ~isa(option, 'char')
        error('option must be a string');
    end
    verbose = strcmp(option, 'verbose');
end

decoder.code = code;
decoder.nbits = nbits;
decoder.type = 'lookup';
decoder.decode = @lookupdecode;

% Build a decode tree for the code.  Useful for fallback when fast decoding is
% not possible:
decoder.slowdecoder = treedecoder(code);
decoder.n_symbols = decoder.slowdecoder.n_symbols;
decoder.table(2 ^ nbits) = struct('isvalid', [], 'indices', [], 'remaining', []);

% Build a table with all possible codemessages with nbits bits:
codemessages = dec2bin(0 : 2 ^ nbits - 1);

hash = ceil((2 ^ nbits) / 79);

% For all possible codemessages, obtain corresponding message and number of 
% remaining bits:
for i = 1 : 2 ^ nbits,
    if verbose
        if mod(i - 1, hash) == 0
            fprintf('#');
        end
    end

% mms 2002/5/22: redundant codes were not possible, not all codemessages are 
% decodeable.  Such codemessages must be flagged as non-decodeable!
% This might be solved using try/catch!  It would be faster (?).
    if iscodemessage(codemessages(i, :), decoder.slowdecoder)
        decoder.table(i).isvalid = true;
        [indices, remaining] = feval(decoder.slowdecoder.decode, ...
            codemessages(i, :), decoder.slowdecoder);
        decoder.table(i).indices = indices;
        decoder.table(i).remaining = length(remaining);
    else
        decoder.table(i).isvalid = false;
    end
end

if verbose
    fprintf('\n');
end
