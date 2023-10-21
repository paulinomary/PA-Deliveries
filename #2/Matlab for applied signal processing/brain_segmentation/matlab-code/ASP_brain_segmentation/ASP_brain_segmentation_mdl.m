% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the Minimum Description Length of the model
% parameters. 
% To increase the speed of the convergence and to simplify the final
% model we divide L by 10

function MDL_out=mdl(VarTheta, VarNbrComp, VarNbrBin, VarSample)

MDL_out = 0;						% MDL value
L = 0;									% ML of the histogram
apriori = log10(2)*64; 	% precision of one double log10

for k=1:VarNbrBin
	tmp = 0;
	for j=1:VarNbrComp
		tmp = tmp + gauss(VarTheta(j,1), VarTheta(j,2),k)*VarTheta(j,3);
	end
	L = L - VarSample(k)*log10(tmp+eps);
end

MDL_out = L/10 + VarNbrComp*(3*apriori);
