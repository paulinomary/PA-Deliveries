% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the maximum likelihood of VarNbrComp Gaussian curves
% for any input image whose samples are stored in the histogram VarSample.
% It uses a smart property of our special case: Instead of computing the 
% probability for each sample (i.e. image pixel), we compute once the
% probability and weight its value with the appearance of the pixel in the
% image. It reduces dramaticaly the computation time (for example 256 
% probabilities are computer in stead of more than 30.000).
% The equation numbers are from the book of Duda: "Pattern Classification
% and Scene Analysis".

function [VarTheta,VarSample]=ASP_brain_segmentation_maxlikelilood(VarTheta, VarNbrComp, VarSample, VarNbrBin, NbrSamples,VarBinSize)

precision = 0.05; % Digit precision of mu and sigma
stop = 0; 	  % Stopping criterion
%VarBinSize = 1;

for i=1:VarNbrComp
	Theta_old(i,1)=VarTheta(i,1); % mu
	Theta_old(i,2)=VarTheta(i,2); % sigma / 10
	Theta_old(i,3)=VarTheta(i,3); % P(wi)
end

while stop==0
	for i=1:VarNbrComp
		Sum_P = 0;
		Sum_P_x = 0;
		Sum_P_x_mu = 0;
		for k=1:VarNbrBin
			P_tmp = 0;
			for j=1:VarNbrComp % Denominator of Eqn (17)
     		   %P_tmp = P_tmp + gauss(k,Theta_old(j,1), Theta_old(j,2))*Theta_old(j,3);
               P_tmp = P_tmp + ASP_brain_segmentation_gauss(k*VarBinSize,Theta_old(j,1), Theta_old(j,2))*Theta_old(j,3);
			end
			if P_tmp == 0
				P_w_x_t = 1;
			else
				%P_w_x_t = VarSample(k)*gauss(k,Theta_old(i,1),Theta_old(i,2))*Theta_old(i,3)/P_tmp;
                P_w_x_t = VarSample(k)*ASP_brain_segmentation_gauss(k*VarBinSize,Theta_old(i,1),Theta_old(i,2))*Theta_old(i,3)/P_tmp;
			end
			Sum_P = Sum_P + P_w_x_t; 	% Eqn 14
			%Sum_P_x = Sum_P_x + P_w_x_t*k;	% Eqn 15
            Sum_P_x = Sum_P_x + P_w_x_t*k*VarBinSize;
			Sum_P_x_mu = Sum_P_x_mu + P_w_x_t*(k*VarBinSize-Theta_old(i,1))*(k*VarBinSize-Theta_old(i,1));
            %Sum_P_x_mu = Sum_P_x_mu + P_w_x_t*(k-Theta_old(i,1))*(k-Theta_old(i,1));% Eqn 16
        end

       VarTheta(i,1) = Sum_P_x/Sum_P;
       VarTheta(i,3) = Sum_P/NbrSamples;
       Sum_P = 0;
	   Sum_P_x = 0;
	   Sum_P_x_mu = 0;
           for k=1:VarNbrBin
			P_tmp = 0;
			for j=1:VarNbrComp % denominator of Eqn (17)
            	%P_tmp = P_tmp + gauss(k,Theta_old(j,1), Theta_old(j,2))*Theta_old(j,3);
                P_tmp = P_tmp + ASP_brain_segmentation_gauss(k*VarBinSize,Theta_old(j,1), Theta_old(j,2))*Theta_old(j,3);
			end
			if P_tmp == 0
				P_w_x_t = 1;
			else
				% P_w_x_t = VarSample(k)*gauss(k,Theta_old(i,1),Theta_old(i,2))*Theta_old(i,3)/P_tmp;
                P_w_x_t = VarSample(k)*ASP_brain_segmentation_gauss(k*VarBinSize,Theta_old(i,1),Theta_old(i,2))*Theta_old(i,3)/P_tmp;
			end
			Sum_P = Sum_P + P_w_x_t; 	% Eqn 14
			% Sum_P_x = Sum_P_x + P_w_x_t*k;	% Eqn 15
			% Sum_P_x_mu = Sum_P_x_mu + P_w_x_t*(k-VarTheta(i,1))*(k-VarTheta(i,1));	% Eqn 16
            Sum_P_x = Sum_P_x + P_w_x_t*k*VarBinSize;
			Sum_P_x_mu = Sum_P_x_mu + P_w_x_t*(k*VarBinSize-Theta_old(i,1))*(k*VarBinSize-Theta_old(i,1));
		end
		 	    
		VarTheta(i,2) = sqrt(Sum_P_x_mu/Sum_P); % Eqn 16
		if VarTheta(i,2) < 1
			VarTheta(i,2) = 1;
        end
        
		 	% Eqn 14
    end
    
	%VarTheta --- stopping criterion?
	for i=1:VarNbrComp
		if abs(VarTheta(i,1)-Theta_old(i,1)) < precision
			stop = stop + 1;
		end
		if abs(VarTheta(i,2)-Theta_old(i,2)) < precision
			stop = stop + 1;
		end
    end
    
	if stop == 2*VarNbrComp
		stop = 1;
	else
		stop = 0;
	end
	% ---
    
	for i=1:VarNbrComp
		Theta_old(i,1)=VarTheta(i,1)
		Theta_old(i,2)=VarTheta(i,2)
		Theta_old(i,3)=VarTheta(i,3)
	end
end
	
