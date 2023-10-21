function imr = dct_q_iq_idct8(im,step)
imr = uint8(idct2(iquantizedct8(quantizedct8(dct2(im),step),step)));
