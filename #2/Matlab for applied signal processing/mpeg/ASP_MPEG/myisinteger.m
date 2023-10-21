function r = myisinteger(m)

% myisinteger - True for matrices containing integers.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

if(~isnumeric(m))
    r = false;
    return
end

if(isa(m, 'double') | isa(m, 'single'))
    frac = m - floor(m);
    r = all(frac(:) == 0 | m(:) == Inf | m(:) == -Inf);
else
    r = true;
end

