function imshow_merror(e,t)
if (nargin<2),
    t = 255*2;
end
e2 = (double(e)/t) + 0.5;
imshow(e2)