% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the value of a Gaussian distribution whitd parameters mu and
% sigma at the point (or vector) x.
% function y=gauss(x,mu,sigma)


function y=ASP_brain_segmentation_gauss(x,mu,sigma)

y = exp(-0.5*(x-mu).*(x-mu)/(sigma*sigma))/(sqrt(2*pi)*sigma);
