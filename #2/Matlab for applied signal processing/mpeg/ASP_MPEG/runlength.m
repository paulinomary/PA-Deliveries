function [nz, vv] = runlength(v, len)

%runlength   Obtain zero run lengths and values
%
%   [nz, vv] = runlength(v) returns the lengths of the runs of zeros in v and 
%   the values of v which are non-zero.  If vector 'v' ends with at least one
%   zero, the last run length is set to -1.
%
%   [nz, vv] = runlength(v, len) returns the concatenation of runlength
%   results when applied to segments of v of length len.
%
%   Example:
%       [nz, vv] = runlength([1 2 0 0 0 3 0 4 5])
%           nz = [0 0 3 1 0]
%           vv = [1 2 3 4 5]
%       [nz, vv] = runlength([1 2 0 0 0 3 0 4 5 0 0 0])
%           nz = [0 0 3 1 0 -1]
%           vv = [1 2 3 4 5]
%
%   See also irunlength, zigzag, izigzag.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if(~isrow(v) | ~isnumeric(v))
    error('v must be a numeric row vector');
end

if nargin >= 2
    if ~isscalar(len) | ~myisinteger(len) | mod(length(v), len) ~= 0
        error('len must be an integer divisor of the length of v');
    end
else
    len = length(v);
end

nz = zeros(1, len);
vv = zeros(1, len);
    
inz = 1;
ivv = 1;

for i = 1 : length(v) / len
    vp = v((i - 1) * len + 1 : i * len);
	vd = double(vp);
	
	% Make sure it terminates with a non-zero.
	vd = [vd 1];
	z = vd == 0;
	cz = cumsum(z);
	nzp = cz(~z);
	nzp = nzp - [0 nzp(1 : end - 1)];
	if nzp(end) == 0
        nzp = nzp(1 : end - 1);
	else
        nzp(end) = -1;
	end
	
	z = z(1 : end - 1);
	vvp = vp(~z);
    
    pinz = inz + length(nzp);
    pivv = ivv + length(vvp);
    nz(inz : pinz - 1) = nzp;
    vv(ivv : pivv - 1) = vvp;
    inz = pinz;
    ivv = pivv;
end

nz(pinz : end) = [];
vv(pivv : end) = [];
