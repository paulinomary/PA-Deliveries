function m = meanpsnr(p)

mse = 255*255 ./ (10.^(p./10));
m = 10*log10(255*255/mean(mse));
