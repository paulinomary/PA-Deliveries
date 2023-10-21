% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function output is a histogram of a 3D MRI image stored in a matrix.
%
% Usage:
%
% Input A is the Matrix containing the 3D image data
% Input loffset is the left offset that is used to remove the background black
% Input roffset is the right offset to cut off the histogram (if roffset=0 no
%  no offset is used.
% Output H is a vector containing the histogram data
% Output X contains the position of the bin centers.
% Output nBin contains the number of bins used
%
% The number of bins will be the same as the number of gray levels minus
% the offsets.
% So the 3D matrix must contain integers (whole numbers).

function [H,X,nBin,binSize] = ASP_brain_segmentation_histogram(A,loffset,roffset)

% get the size of the 3D data
s=size(A);

% Maximum intensity level in the volume
n=floor(max(max(max(A)))-1);

%Define Bins of the image histogram to be used
step=1;
if roffset==0
	X=loffset:step:n;
else
	X=loffset:step:(n-roffset);
end

% ***--- create the historgram ---***

% then calculate the histogram
H=histc(A(:),X);

% the number of bins is the size of X
nBin=size(X);
nBin=nBin(1,2);

binSize=step;

% ***--- end ---***
