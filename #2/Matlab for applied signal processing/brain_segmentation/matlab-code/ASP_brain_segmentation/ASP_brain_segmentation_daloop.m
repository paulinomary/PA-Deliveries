% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performs de MAP using only one iteration of ICM

function [nL]=ASP_brain_segmentation_daloop(L,Image,mu,sigma,NbrComp,Beta);

s=size(Image);
T=zeros(NbrComp,s(1),s(2),s(3));

% Voxel size dimensions
Sx=1.;
Sy=1.;
Sz=1.3;

for l=1:NbrComp,
    disp('Computing energy UX + UYX...');
    % Reconstruction of 4D matrix: every element contains a 5 element vector representing each label energy.
    if NbrComp==5,
        %T(l,:,:,:)=calcuyx(mu(l),sigma(l),Image)+0.6*calcux_5class(L,Sx,Sy,Sz,l);
        T(l,:,:,:)=ASP_brain_segmentation_calcuyx(mu(l),sigma(l),Image)+Beta*ASP_brain_segmentation_calcux_5class(L,Sx,Sy,Sz,l);
    elseif NbrComp==3,
        T(l,:,:,:)=ASP_brain_segmentation_calcuyx(mu(l),sigma(l),Image)+Beta*ASP_brain_segmentation_calcux_3class(L,Sx,Sy,Sz,l); 
    end
    clear UX
    clear UYX
end

% Min value position corresponds to the new label
[V,newLabel]=min(T);
clear T
nL=reshape(newLabel,s(1),s(2),s(3));
nL(find(L==0))=0;
