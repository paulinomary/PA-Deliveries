function r = prod(c)

% for now works only for cell vectors

r = 1;

for i = 1 : length(c),
    r = r .* c{i};
end