% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates P(x) given a label l and a location
% in the labelmap L x,y,z

function [PlX]=ASP_brain_segmentation_calcpx(L,Sx,Sy,Sz,label,NbrComp) 

%Normalization factor (maximum value)
Ux_max=2/Sx+2/Sy+2/Sz+4/sqrt(Sx.^2+Sy.^2);
Norm=exp(0.6*Ux_max);

if (NbrComp==5),
    Ux=ASP_brain_segmentation_calcux_5class(L,Sx,Sy,Sz,label);
    PlX=(1/Norm).*(exp(-0.6*Ux));
elseif (NbrComp==3),
    Ux=ASP_brain_segmentation_calcux_3class(L,Sx,Sy,Sz,label);
    PlX=(1/Norm).*(exp(-0.6.*Ux));
end
clear Ux

s=size(PlX);
NL=s(1);
NC=s(2);
NS=s(3);

PlX(:,:,1)=ones(NL,NC);
PlX(:,:,NS)=ones(NL,NC);
PlX(1,:,:)=ones(NC,NS);
PlX(NL,:,:)=ones(NC,NS);
PlX(:,1,:)=ones(NL,NS);
PlX(:,NC,:)=ones(NL,NS);

PlX(find(L==0))=0;


