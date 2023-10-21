% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to calculate Pt(l|yi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Where Image is the 3D brain image
% L is the label map
% x, y and z are the coordinates of the point
% l is the label that is considered
% nH is the normalized histogram
% X are the histogram bin centers
% meanl is the estimated mean for label l
% sigmal is the estimated sigma for label l


function [ptly]=ASP_brain_segmentation_calcptly(Image, L, l, nH, X, meanl, sigmal,NbrComp);

Sx=1.;
Sy=1.;
Sz=1.3;

disp('Calculating P(X)...');
PlX=ASP_brain_segmentation_calcpx(L,Sx,Sy,Sz,l,NbrComp);

disp('Calculating P(Y|X)...');
Gt=ASP_brain_segmentation_gauss(Image,meanl,sigmal);

disp('Calculating P(Y)...');
pyi=nH(round(Image-X(1))); % chance of gray level yi is the histogram value

ptly=zeros(size(PlX));
ptly=(Gt.*PlX)./pyi;

% to overcome problems when pyi=0, don't know if this is the solution
ptly(find(pyi==0))=0;

% to not consider intensity points under 2
ptly(find(Image<2))=0;

% if l==1,
%     save pCSF ptly;
% else if l==2,
%         save pGM ptly;
%     else
%         save pWM ptly;
%     end
% end

clear PlX
clear Gt
clear pyi



