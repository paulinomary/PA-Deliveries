% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.

% 3D Gaussian filtering (separable implementation)
function [sV]=ASP_brain_segmentation_smooth_3DGauss(filename,Nx,Ny,Nz,dataType,SizeFilter,Sigma);

% Open the image to smooth
fid=fopen(filename,'r','l');
V=fread(fid,Nx*Ny*Nz,dataType');
V=reshape(V,Nx,Ny,Nz);

% Smooth data with a 3D Gaussian filter
sV=smooth3(V,'gaussian',SizeFilter,Sigma);
