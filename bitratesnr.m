

x1=84.51;
x2=38.85;
x3=26.36;
x4=13.63;



y1=8000*16;
y2=8000*8;
y3=8000*6;
y4=8000*4;


figure
plot(x1,y1,'x',x2,y2,'x',x3,y3,'x',x4,y4,'x','LineWidth', 4)
title('Gráfico SNR vs BitRates para o áudio outfemale')
ylabel('bitrate') 
xlabel('SNR') 
legend('16 bits', '8 bits', '6 bits', '4 bits','Location','northwest','Orientation','vertical')
