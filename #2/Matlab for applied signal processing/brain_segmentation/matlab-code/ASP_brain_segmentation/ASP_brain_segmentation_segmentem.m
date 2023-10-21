% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function will segment an image using an EM algorithm to 
% fit a number Gaussians on the histogram of a 3D MRI data set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image is the input 3D MRI data
% loffset is the left offset to be used to calculate the histogram
% 	(i.e. use it to cut off the large peak at low gray levels in the
% 	image background)
% roffset is the right offset
% Theta is a matrix containing the estimates of mu sigma and p(wi)
% 	The form is like this:
%		[	mu1 sigma1 pw1
%			mu2 sigma2 pw2
%			mu3 sigma3 pw3]
%
% S is the segmented image
% e are the used boundaries
% Theta is the newly found Theta
% G are the found gaussians
% fitper is the percentage of fit of the histogram

function [S,e,Theta,fitper,G]=ASP_brain_segmentation_segmentem(Image, loffset, roffset, Theta, binSize)

s=size(Theta);
NbrComp=s(1);
clear s;

% Computing normalized image histogram
[H,X,NbrBin,BinSize]=ASP_brain_segmentation_histogram(Image,loffset,roffset);
NbrSamples=sum(sum(H));
nH=H*1/NbrSamples;

stop = 0;

% just initilizing
MDL_old=2147483648;
MDL_new=2147383648;
ThetaOld=zeros(size(Theta));
n=0;

% Find out the new Theta
while stop==0
  if MDL_old<=MDL_new
    stop = 1;
  else
    disp('Maxlikelihood');
    [Theta, Samplem]=ASP_brain_segmentation_maxlikelihood(Theta, NbrComp, H, NbrBin,NbrSamples, binSize);
    ThetaOld=Theta;
    MDL_old=MDL_new;
    MDL_new=ASP_brain_segmentation_mdl(Theta,NbrComp,NbrBin,H);
    disp('Difference of Minimum Descriptor Lenght')
    MDL_old-MDL_new
		
		% --- extra stopping criterion because the other won't stop.
        disp('Theta-ThetaOld')
		t1=sum(sum(abs(Theta-ThetaOld)))
		
		if t1==0;
            disp('Iteration')
			n=n+1
		end
		if n>50;
			stop =1;
		end
		% --- end of the criterion
  end
end

for i=1:NbrComp
	Theta(i,1)=((loffset-1)+Theta(i,1));
	Theta(i,2)=Theta(i,2);
end

Theta=sortrows(Theta); % sort from small to large mu.

% find the boundaries of gray levels (crossings of the Gaussians)
% calculate the Gaussians
for i=1:NbrComp
	for j=1:length(X)
  		G(i,j)=binSize*Theta(i,3).*ASP_brain_segmentation_gauss(X(j),Theta(i,1),Theta(i,2));
	end
end

% find the boundaries of gray levels

%debug
sg=size(G);
sg=sg(1);

e(1)=0;

if sg==5
	g=double(G(2,:)>G(1,:));
	g=g+double(G(3,:)>G(2,:));
	g=g+double(G(4,:)>G(3,:));
	g=g+double(G(5,:)>G(4,:));
elseif sg==4
	g=double(G(2,:)>G(1,:));
        g=g+double(G(3,:)>G(2,:));
        g=g+double(G(4,:)>G(3,:));
elseif sg==3
        g=double(G(2,:)>G(1,:));
        g=g+double(G(3,:)>G(2,:));
end

l=2;
for i=1:length(X)-1
        if g(i+1)>g(i)
                e(l)=X(i); % gray level boundaries
                l=l+1;
        end
end

clear s;
s=size(e);
hmax=max(max(nH));

% plot the fitted gaussians and boundaries
figure;
plot(X,nH,':');
hold on;
for i=1:NbrComp
  plot(X,G(i,:),'--');
end

plot(X,sum(G),'-');
hold on;

for i=1:s(2)
	line([e(i) e(i)],[0 hmax+hmax/10]);
end

title('Histogram with the fitted Gaussians and boundaries')
xlabel('Gray level')

% squared difference between histogram and fit.

fitper=100*(NbrSamples-sum(abs(sum(G)*NbrSamples-H')))/(NbrSamples)

pause(0.1)

% doing the actual segmentation
S=zeros(size(Image));

for i=1:s(2);
  S=S+double(Image>e(i));
end

