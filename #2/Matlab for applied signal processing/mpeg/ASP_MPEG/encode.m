function codemessage = encode(message, alphabet, code)

%encode   Message encoding.
%
%   codemessage = encode(message, alphabet, code) returns the encoded version of 
%   message 'message' consisting of symbols from the alphabet 'alphabet'.  The 
%   code used is given by 'code'.
%
%   'message' and 'alphabet' are row cell vectors of symbols.  Any class of 
%   symbol is acceptable provided it has the method 'equal' overloaded.
%
%   codemessage = encode(indexmessage, code) returns the encoded version of 
%   index message 'indexmessage' consisting of integer symbols in the range 
%   1 : length(code).  I.e., the index message elements are indices into an
%   implicit alphabet.  The code used is given by 'code'.
%
%   'indexmessage' is a row vector of integers in the range 1 : length(code).
%
%   In both cases code 'code' must be a cell of row char vectors (strings).  Each
%   element of this array is the binary codeword for the corresponding symbol in 
%   the alphabet.
%
%   The code message returned is a string with the concatenation of the codewords
%   corresponding to the message symbols.
% 
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if(nargin ~= 2 & nargin ~= 3)
    error('must have 2 or 3 arguments: message, alphabet and codes or indexmessage and code');
end

if(nargin == 3)
    if(~iscell(message) | ~iscell(alphabet) | ~iscell(code))
        error('message, alphabet and code must all be cell arrays');
    elseif(~isrow(message) | ~isrow(alphabet) | ~isrow(code))
        error('message, alphabet and code must all be row cell vectors');
    elseif(length(alphabet) ~= length(code))
        error('alphabet and codes must have same length');
    end
else
    indexmessage = message;
    code = alphabet;
    if(~iscell(code))
        error('code must all a cell array');
    elseif(~myisinteger(indexmessage))
        error('index message must be an integer array');
    elseif(~isrow(indexmessage) | ~isrow(code))
        error('index message and code must both be row cell vectors');
    end
end
    
if(nargin == 3)
    indexmessage = messagetoindex(message, alphabet);
end

% Simple, but slow (for long messages! for short messages this version is
% faster):
%codemessage = char([code{indexmessage}]);

% Build an array of bits:
bits = char(code(indexmessage))';

% Reshape it to obtain the coded message, though still with spaces for padding:
codemessage = bits(:)';

% Remove padding:
codemessage(codemessage == ' ') = [];
