% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates posterior probability maps for each tissue
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Where Image is the 3D brain image
% L is the label map
% mu is a vector containing the estimated means for each label 
% sigma is a vector containing the estimated sigma's for each label
% NbrComp is the number of classes

function []=ASP_brain_segmentation_ProbMaps(Image, L, mu, sigma, NbrComp, FileName);

% IMAGE HISTOGRAM
% ---------------
Loffset=1; % Comment this line if Loffset is chosen by the user
Roffset=0; % Comment this line if Roffset is chosen by the user
[H,X,NbrBin]=ASP_brain_segmentation_histogram(Image,Loffset,Roffset);
NbrSamples=sum(sum(H));
% Normalized histogram
nH=H*1/NbrSamples;
% ---------------

s=size(Image);
Backup=Image; % Backup of original Image
index=find(Image<(X(1)+1)); % To avoid (Image-X(1))=0 in calcptly.m
Image(index)=(X(1)+1);
ptly=zeros(NbrComp,s(1),s(2),s(3));

for l=1:NbrComp
    disp('Calculating P(Y|X)...');
    ptly(l,:,:,:)=ASP_brain_segmentation_calcptly(Image, L, l, nH, X, mu(l), sigma(l), NbrComp); % Equation 27 (Y. Zhang article)
end

Image=Backup;
clear Backup

%Normalization of the probability maps
disp('Normalization of the probability maps (the sum of all tissue probabilities must be equal 1)');

warning off MATLAB:divideByZero
a=ptly(1,:,:,:);
A=reshape(a,s(1),s(2),s(3));
A=A./max(max(max(A)));
a=ptly(2,:,:,:);
B=reshape(a,s(1),s(2),s(3));
B=B./max(max(max(B)));
a=ptly(3,:,:,:);
C=reshape(a,s(1),s(2),s(3));
clear a
C=C./max(max(max(C)));

D=zeros(s(1),s(2),s(3));
E=D;

if NbrComp==5,
    a=ptly(4,:,:,:);
    D=reshape(a,s(1),s(2),s(3));
    D=D./max(max(max(B)));
    a=ptly(5,:,:,:);
    E=reshape(a,s(1),s(2),s(3));
    clear a;
    E=E./max(max(max(C)));
end
clear ptly;

if NbrComp==3,
    Total=A+B+C;
    A=A./Total;
    B=B./Total;
    C=C./Total;
    A(find(Total==0))=0;
    B(find(Total==0))=0;
    C(find(Total==0))=0;
    
    %write output files
    fname=strcat(FileName,'_Pcsf.raw');
    fid = fopen(fname,'wb');
    fwrite(fid,A*255,'float');
    fclose(fid);
    fname=strcat(FileName,'_Pgm.raw');
    fid = fopen(fname,'wb');
    fwrite(fid,B*255,'float');
    fclose(fid);
    fname=strcat(FileName,'_Pwm.raw');
    fid = fopen(fname,'wb');
    fwrite(fid,C*255,'float');
    fclose(fid);
    
    figure; imagesc(squeeze(A(:,99,:)));title('Probability of CSF'); 
    figure; imagesc(squeeze(B(:,99,:)));title('Probability of GM'); 
    figure;  imagesc(squeeze(C(:,99,:)));title('Probability of WM');
    
elseif NbrComp==5,
    Total=A+B+C+D+E;
    A=A./Total;
    B=B./Total;
    C=C./Total;
    D=D./Total;
    E=E./Total;
    A(find(Total==0))=0;
    B(find(Total==0))=0;
    C(find(Total==0))=0;
    D(find(Total==0))=0;
    E(find(Total==0))=0;
    
    %write output files
    fname=strcat(FileName,'_Pcsf.raw');
    fid = fopen(fname,'wb');
    fwrite(fid,A*255,'float');
    fclose(fid);
    fname=strcat(FileName,'_Pcg.raw');
    fid = fopen(fname,'wb');
    fwrite(fid,B*255,'float');
    fclose(fid);
    fname=strcat(FileName,'_Pgm.raw');
    fid = fopen(fname,'wb');
    fwrite(fid,C*255,'float');
    fclose(fid);
    fname=strcat(FileName,'_Pgw.raw');
    fid = fopen(fname,'wb');
    fwrite(fid,D*255,'float');
    fclose(fid);
    fname=strcat(FileName,'_Pwm.raw');
    fid = fopen(fname,'wb');
    fwrite(fid,E*255,'float');
    fclose(fid);
end

warning on MATLAB:divideByZero
clear A
clear B
clear C
clear D
clear E


