% Picanha

%% 1. Usando a função imread do MATLAB leia as imagens a usar (ver nota no fim) para memória e apresente-as no écran.

picanha = imread("picanha.png");
%image(picanha);
figure, imshow(picanha), title('Original Image')

picanhaYCbCr = rgb2ycbcr(picanha);

%% 2. Usando a função MATLAB rgb2ycbcr converta as imagens de RGB (formato original) para YCbCr. Exibir as componentes R, G e B da imagem original
subplot(2,3,1);
imshow(picanha(:,:,1));
title('Componente R');

subplot(2,3,2);
imshow(picanha(:,:,2));
title('Componente G');

subplot(2,3,3);
imshow(picanha(:,:,3));
title('Componente B');

% Componentes Y, Cb e Cr
subplot(2,3,4);
imshow(picanhaYCbCr(:,:,1));
title('Componente Y (YCbCr)');

subplot(2,3,5);
imshow(picanhaYCbCr(:,:,2));
title('Componente Cb (YCbCr)');

subplot(2,3,6);
imshow(picanhaYCbCr(:,:,3));
title('Componente Cr (YCbCr)');

%% 3. Pretende-se agora verificar empiricamente a importância subjectiva relativa das componentes Y, Cb e Cr forçando uma quantização por eliminação da informação menos significativa de cada pixel.
% a) b) c)

Y = picanhaYCbCr(:,:,1);
new_Y = 8*floor(Y/8);
new_img= picanhaYCbCr;
new_img(:,:,1) = new_Y;
%image(new_img);
figure, imshow(new_img), title('YCbCr only Y 3 bits')

RGB1 = ycbcr2rgb(new_img);
%image(RGB1);
figure, imshow(RGB1), title('RGB YCbCr only Y 3 bits')

Cb = picanhaYCbCr(:,:,2);
new_Cb = 8*floor(Cb/8);
new_image2 = picanhaYCbCr;
new_image2(:,:,2) = new_Cb;
%image(new_image2);
figure, imshow(new_image2), title('YCbCr only Cb 3 bits')


RGB2 = ycbcr2rgb(new_image2);
%image(RGB2);
figure, imshow(RGB2), title('RGB YCbCr only Cb 3 bits')


Cr = picanhaYCbCr(:,:,3);
new_Cr = 8*floor(Cr/8);
new_image3 = picanhaYCbCr;
new_image3(:,:,3) = new_Cr;
%image(new_image3);
figure, imshow(new_image3), title('YCbCr only Cr 3 bits')

RGB3 = ycbcr2rgb(new_image3);
%image(RGB3);
figure, imshow(RGB3), title('RGB YCbCr only Cr 3 bits')

%% 3. D)Repita a) a c) colocando a zero 4 bits

new_Y2 = 16*floor(Y/16);
new_img5= picanhaYCbCr;
new_img5(:,:,1) = new_Y2;
%image(new_img);
figure, imshow(new_img5), title('YCbCr only Y 4 bits')

RGB5 = ycbcr2rgb(new_img5);
%image(RGB3);
figure, imshow(RGB5), title('RGB YCbCr only Y 4 bits')

new_Cb2 = 16*floor(Cb/16);
new_img6= picanhaYCbCr;
new_img6(:,:,1) = new_Cb2;
%image(new_img);
figure, imshow(new_img6), title('YCbCr only Cb 4 bits')

RGB6 = ycbcr2rgb(new_img6);
%image(RGB3);
figure, imshow(RGB6), title('RGB YCbCr only Cb 4 bits')

new_Cr2 = 16*floor(Cr/16);
new_img7= picanhaYCbCr;
new_img7(:,:,1) = new_Cr2;
%image(new_img);
figure, imshow(new_img7), title('YCbCr only Cr 4 bits')

RGB7 = ycbcr2rgb(new_img7);
%image(RGB3);
figure, imshow(RGB7), title('RGB YCbCr only Cr 4 bits')


%%
% 5. Gere de novo as componentes Y das imagens de teste (partindo sempre das imagens originais). Usando a função imnoise, adicione ruído gaussiano de média nula e desvio 0.025 à componente Y.

new_Y2 = imnoise(Y,'gaussian',0,0.025);
new_img4= picanhaYCbCr;
new_img4(:,:,1) = new_Y2;
%image(new_img4);
figure, imshow(new_img4), title('Y with gaussian noise')

RGB4 = ycbcr2rgb(new_img4);
%image(RGB4);
figure, imshow(RGB4), title('RGB with gaussian noise')


