function [message, codemessage] = decodemessage(codemessage, alphabet, code, n, option)

%decodemessage   Message decoding.
%
%   [message, codedmessage] = decodemessage(codemessage, alphabet, code[, n[, option]]) 
%   returns the decoded version of the coded message 'codemessage' consisting of
%   a sequence of codewords from the prefix code 'code'.  The symbols are taken
%   from alphabet 'alphabet'.
%
%   'alphabet' is a row cell vectors of symbols.  Any class of symbol is 
%   acceptable provided it has the method 'equal' overloaded.
%
%   The returned message 'message' corresponds to a row cell of symbols from the
%   alphabet.
%
%   [indexmessage, codemessage] = decodemessage(codemessage, code[, n[, option]]) 
%   returns the decoded version of the coded message 'codemessage' consisting of
%   a sequence of codewords from the prefix code 'code'.  The decoded message is 
%   returned in the form of indices into an implicit alphabet. 
%
%   The returned index message 'indexmessage' corresponds to a row cell of 
%   indices into an implicit alphabet.
%
%   Instead of 'code', a decoder may be passed, such as returned by functions
%   'lookupdecoder' and 'treedecoder'.  It is recomended to calculate the decoder
%   beforehand so that it may be reused whenever necessary.
%
%   In all cases 'codemessage' is a row char vector (i.e., a string) consisting
%   of '0' and '1'.
%
%   In both cases code 'code' must be a cell of char (row) arrays.  Each element
%   of this array is the binary codeword for the corresponding symbol in the
%   alphabet.
%
%   If given, n is the number of symbols to decode.
%
%   Returns also the undecoded bits of the message, so that they can be decoded
%   later.
%
%   New types of decoders may be created.  Two functions are necessary.  The
%   first creates a decoder from a code (as 'treedecoder').  The second 
%   decodes an index message from a coded message, the created decoder, and
%   optionaly the number of symbols to decode (as 'treedecode').  The decoder
%   must be a scalar struct with the following fields:
%       a) type - a string naming the type of decoder.
%       b) n_symbols - the number of symbols with codewords in the alphabet.  
%          I.e., the number of symbols with non-zero probability.
%       c) decode - a function handle to the specialized decoding function, as
%          decribed.
%       d) Any other field necessary for the decoder.
%
%   Hashmarks are shown during decoding if the option 'verbose' is passed.
%
%   See also lookupdecode, treedecode, lookupdecoder, treedecoder.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if nargin ~= 2 & nargin ~= 3 & nargin ~= 4 & nargin ~= 5
    error(['must have 2 to 5 arguments: codedmessage, alphabet, ' ...
        'codes/decoder [, n] or codedmessage, codes/decoder[, n[, option]]']);
end

if(nargin >= 3 & (iscell(code) | isstruct(code)))
    if(nargin ~= 3 & nargin ~= 4  & nargin ~= 5)
        error('wrong number of arguments with a given alphabet');
    end
    if(nargin >= 4)
        if(~isscalar(n) | ~isinteger(n))
            error('n must be a scalar integer');
        elseif(n < 0)
            error('n must be non-negative');
        end
        ngiven = true;
    else
        ngiven = false;
    end
    
    if nargin < 5
        verbose = false;
    else
        if ~isa(option, 'char')
            error('option must be a string');
        end
        verbose = strcmp(option, 'verbose');
    end

    if(~iscell(alphabet))
        error('alphabet and code must both be cell arrays');
    elseif(~isbinarychar(codemessage))
        error('codemessage must be a binary string');
    elseif(~isrow(codemessage) | ~isrow(alphabet))
        error('codemessage and alphabet must both be row cell vectors');
    end
    if(iscell(code))
        if(~isrow(code))
            error('code must be a row cell vector');
        elseif(length(alphabet) ~= length(code))
            error('alphabet and codes must have same length');
        end
    end
else
    if nargin ~= 2 & nargin ~= 3 & nargin ~= 4
        error('wrong number of arguments without an alphabet');
    end

    if nargin < 4
        verbose = false;
    else
        option = n;
        if ~isa(option, 'char')
            error('option must be a string');
        end
        verbose = strcmp(option, 'verbose');
    end

    if nargin >= 3
        n = code;
        if(~isscalar(n) | ~myisinteger(n))
            error('n must be a scalar integer');
        elseif(n < 0)
            error('n must be non-negative');
        end
        ngiven = true;
    else
        ngiven = false;
    end

    code = alphabet;
    alphabet = cell(1, 0);
    if(~iscell(code) & ~isstruct(code))
        error('code must be either be a cell array or a decoder struct');
    elseif(iscell(code))
        if(~isrow(code))
            error('code must be a row cell vector');
        end
    end
    if(~isbinarychar(codemessage))
        error('codemessage must be a binary string');
    elseif(~isrow(codemessage) | ~isrow(code))
        error('codemessage and code must both be row vectors');
    end
end
    
if verbose
    verboseargs = {'verbose'};
else
    verboseargs = {};
end

% Build a decoder, if necessary.  Uses table lookup by default:
    
if(isstruct(code))
    decoder = code;
    clear code;
else
    if length(codemessage) < 2^10*10
        if verbose
            fprintf('Generating tree decoder:\n');
        end
        decoder = treedecoder(code, verboseargs{:});
    else
        if verbose
            fprintf('Generating lookup decoder with 10 bits:\n');
        end
        decoder = lookupdecoder(code, 10, verboseargs{:});
    end
end
   
%fprintf('using a %s decoder: ', decoder.type);
%t0 = cputime;

% Decode:

%  Deal with special case where there are no symbols!
if decoder.n_symbols == 0
	if isempty(alphabet)
        message = zeros(1, 0);
	else
        message = cell(1, 0);
	end
    return;
end

%  Deal with special case where there is a single symbol with probability 1.  
%  In that case the number of symbols to decode 'n' must be given!
if decoder.n_symbols == 1
    if ~ngiven
        error('cannot guess number of symbols in message if alphabet is degenerate');
    end
	if isempty(alphabet)
        message = ones(1, n);
	else
        message = cell(1, n);
        message(1:n) = alphabet(1);
	end
    return
end

if ~ngiven
    n = Inf;
end

if verbose
    fprintf('Decoding:\n');
end

[indexmessage, codemessage] = feval(decoder.decode, codemessage, decoder, n,...
    verboseargs{:});

if isempty(alphabet)
    message = indexmessage;
else
    message = indextomessage(indexmessage, alphabet);
end

%fprintf('%f seconds\n', cputime - t0);