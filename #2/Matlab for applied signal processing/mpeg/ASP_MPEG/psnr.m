function e = psnr(i1, i2)

if(isa(i1, 'uint8'))
    max = 255;
else
    max = 1;
end

e = 20 * log10(max / sqrt(sum((double(i1(:)) - double(i2(:))) .^ 2) / prod(size(i1))));