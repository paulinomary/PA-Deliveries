function [e, diff] = mypsnr(i1, i2)

if(isa(i1, 'uint8'))
    max = 255;
else
    max = 1;
end
diff = sum((double(i1(:)) - double(i2(:))) .^ 2);
e = 20 * log10(max / sqrt(diff / numel(i1)));