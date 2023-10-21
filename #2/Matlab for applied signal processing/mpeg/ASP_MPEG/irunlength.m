function [v, nz, vv] = irunlength(nz, vv, len, totallen)

%irunlength   Obtain vector from zero run lengths and values
%
%   [v, nz, vv] = irunlength(nz, vv, len) returns the vector corresponding to
%   the lengths of the runs of zeros nz and to the non-zero values vv.  The 
%   returned vector will have length len.  The remaining zero run lengths and
%   non-zero values are returned.
%
%   [v, nz, vv] = irunlength(nz, vv, len, totallen) returns the vector
%   corresponding to the lengths of the runs of zeros nz and to the
%   non-zero values vv.  The returned vector will have length totallen,
%   but will be "decoded" in segments of len. The remaining zero run
%   lengths and non-zero values are returned.
%
%   Example:
%       [v, nz, vv] = irunlength([0 0 3 1 0], [1 2 3 4 5], 1)
%           v = 1
%           nz = [0 3 1 0]
%           vv = [2 3 4 5]
%       [v, nz, vv] = irunlength([0 0 3 1 0], [1 2 3 4 5], 2)
%           v = [1 2]
%           nz = [3 1 0]
%           vv = [3 4 5]
% 
%       [v, nz, vv] = irunlength([0 0 3 1 0], [1 2 3 4 5], 3)
%           ??? Error using ==> irunlength
%           nz incompatible with len (end missed)
% 
%       [v, nz, vv] = irunlength([0 0 3 1 0], [1 2 3 4 5], 4)
%           ??? Error using ==> irunlength
%           nz incompatible with len (end missed)
% 
%       [v, nz, vv] = irunlength([0 0 3 1 0], [1 2 3 4 5], 5)
%           ??? Error using ==> irunlength
%           nz incompatible with len (end missed)
% 
%       [v, nz, vv] = irunlength([0 0 3 1 0], [1 2 3 4 5], 6)
%           v = [1 2 0 0 0 3]
%           nz = [1 0]
%           vv = [4 5]
% 
%       [v, nz, vv] = irunlength([0 0 3 1 0], [1 2 3 4 5], 7)
%           ??? Error using ==> irunlength
%           nz incompatible with len (end missed)
% 
%       [v, nz, vv] = irunlength([0 0 3 1 0], [1 2 3 4 5], 8)
%           v = [1 2 0 0 0 3 0 4]
%           nz = 0
%           vv = 5
% 
%       [v, nz, vv] = irunlength([0 0 3 1 0], [1 2 3 4 5], 9)
%           v = [1 2 0 0 0 3 0 4 5]
%           nz = Empty matrix: 1-by-0
%           vv = Empty matrix: 1-by-0
% 
%       [v, nz, vv] = irunlength([0 0 3 1 0], [1 2 3 4 5], 1)
%           v = 1
%           nz = [0 3 1 0]
%           vv = [2 3 4 5]
% 
%   To do:
%       Should return something sensible when the length does not correspond 
%       to an integer number of runs.
%
%   See also irunlength, zigzag, izigzag.
%
%Copyright (C) 2001-2003 - Manuel Menezes de Sequeira (ISCTE).

if(~isrow(nz) | ~isa(nz, 'double') | any(nz < -1))
    error('nz must be a double row vector with values >= -1');
elseif(~isrow(vv) | ~isnumeric(vv))
    error('vv must be a numeric row vector');
elseif(~isscalar(len) | ~isa(len, 'double') | len < 0)
    error('len must be a double non-negative scalar');
end

if nargin >= 4
    if ~isscalar(totallen) | ~myisinteger(totallen) | totallen < 0 | mod(totallen, len) ~= 0
            error('len must be an integer divisor of totallen');    
    end
else
    totallen = len;
end

if len == 0
    v = zeros(1, 0);
    nz = nz;
    vv = vv;
    return
end

v = zeros(1, 0);

for i = 1 : totallen / len
    [v((i - 1) * len + 1 : i * len), nz, vv] = irl(nz, vv, len);
end

function [v, rnz, rvv] = irl(nz, vv, len)

nonz = cumsum(nz + 1);

i = 1;
while i <= length(nonz) & nonz(i) < len & nz(i) ~= -1
    i = i + 1;
end

if i > length(nonz)
    error('nz incompatible with len (too short)');
end

if nonz(i) > len
    error('nz incompatible with len (end missed)');
end

rnz = nz(i + 1 : end);
nonz = nonz(1 : i);

if nz(i) == -1
    rvv = vv(i : end);
    nonz(end) = len;
    vv = [vv(1 : i - 1) 0];
else
    rvv = vv(i + 1 : end);
    nz(end) = len;
    vv = vv(1 : i);
end

v(nonz) = vv(1 : i);