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
% This function calculates the new mean and sigma's for the tissues in the hmm em algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Where Image is the 3D brain image
% L is the label map
% nH is the normalized histogram
% X are the histogram bin centers
% U is a vector containing the estimated means for each label 
% S is a vector containing the estimated sigma's for each label

function [ul, sl]=ASP_brain_segmentation_calculsl(Image, L, nH, X, U, S, NbrComp);

s=size(Image);

topul=zeros(1,NbrComp);
bottomul=zeros(1,NbrComp);
topsl=zeros(1,NbrComp);
bottomsl=zeros(1,NbrComp);

Backup=Image; % Backup of original Image
index=find(Image<2); % To avoid index values becoome 0 (no supported in matlab)
Image(index)=2;

for l=1:NbrComp
    disp('Calculating P(Y|label)...');
    ptly=ASP_brain_segmentation_calcptly(Image,L,l, nH, X, U(l), S(l),NbrComp);
        
    topul(l)=sum(sum(sum(ptly.*Image)));              % Equation 7.22: numerador (Bram's report)
    topsl(l)=sum(sum(sum(ptly.*((Image-U(l)).^2))));  % Equation 7.23: numerador (Bram's report)
    bottomul(l)=sum(sum(sum(ptly)));                  % Equation 7.22: denominador (Bram's report)
    bottomsl(l)=bottomul(l);                          % Equation 7.23: denominador (Bram's report)

    clear ptly 
end

ul=topul./bottomul;
sl=(topsl./bottomsl).^0.5;
Image=Backup;
clear Backup
