function next = next(current, digits)

if nargin < 2
    digits = 10;
end

if isscalar(digits)
    digits = digits * ones(size(current));
end

next = current;

i = 1 : length(current);
largest0 = max(i(next ~= digits - 1));
if isempty(largest0)
    next(1 : end) = 0;
else
    next(largest0) = next(largest0) + 1;
    next((largest0 + 1) : end) = 0;
end

