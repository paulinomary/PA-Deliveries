% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.

close all
clear all

% Load image
disp('')
disp('Please select a 3D brain image file')
[bfilename, bpathname] = uigetfile('*.*', 'Open a .mat file');
filenameMRI=strcat(bpathname,bfilename);
load(filenameMRI);

% Computing image histogram
loffset=1; % Avoid the pic of the histogram at 0 gray-level due to background
roffset=0; 
[H,X,NbrBin,BinSize]=ASP_brain_segmentation_histogram(Axial,loffset,roffset);
NbrSamples=sum(sum(H));
nH=H*1/NbrSamples;

% Gauss fitting
mus = [115 345 410]';
s=size(mus);
NbrComp=s(1);
for i=1:NbrComp
	sig(i,1)=5.0; % Small value different from zero
    pw(i,1)=1/NbrComp; % Initially, the assumption that all classes are equiprobable is done
end

Theta=cat(2,mus,sig,pw)
% Starting calculations
disp(' ');
disp('Now starting EM algorithm...');
[S,e,Theta,fitper,G]=ASP_brain_segmentation_segmentem(Axial, loffset, roffset, Theta, BinSize);
 
disp('Estimated means and variance:');
ul=Theta(:,1)
sl=Theta(:,2)

disp('Estimated thresholds for Bayesian decision:');
e
se=size(e);

% Saving Bayes Classification
filenameMRI=filenameMRI;
fname=strcat(filenameMRI(1:end-4),'_S.raw');
fid = fopen(fname,'wb');
fwrite(fid,S,'ushort');
fclose(fid);

% Reading Bayesian classification
fid=fopen('Image1_voi_S.raw','r','l');
%fid=fopen('Image2_voi_S.raw','r','l');
S=fread(fid,186*186*121,'ushort');
S=reshape(S,186,186,121);

disp('Now going to do the HMM-EM');
nL=S;
clear S;
loffset=1;
roffset=0;
[H,X,NbrBin,BinSize]=ASP_brain_segmentation_histogram(Axial,loffset,roffset);
NbrSamples=sum(sum(H));
nH=H*1/NbrSamples;

% HMRF
% nit=input('How many iterations must be done for the HMM: ');
nit=60;
beta=10;
nLold=nL;
Change=zeros(nit,1);
for ii=1:nit
	disp('Now doing the MAP part');
	nL=ASP_brain_segmentation_daloop(nL,Axial,ul,sl,NbrComp,beta);
    disp('New map nL calculated');
    disp('Computing number of voxels that changed');
    Diff=(nLold-nL);
    Change(ii)=length(find(Diff))./length(find(nL))*100;
    nLold=nL;
    disp('Recompute mean and variance');
    [ul,sl]=ASP_brain_segmentation_calculsl(Axial, nL, nH, X, ul, sl,NbrComp);
end


% Output MAP Classification
filenameMRI=filenameMRI;
fname=strcat(filenameMRI(1:end-4),'_nL.raw');
fid = fopen(fname,'wb');
fwrite(fid,nL,'ushort');
fclose(fid);

% Visualization
% figure
% imagesc(squeeze(nL01(:,54,:)))
% title('Beta 0.1')
% figure
% imagesc(squeeze(nL06(:,54,:)))
% title('Beta 0.6')
% figure
% imagesc(squeeze(nL12(:,54,:)))
% title('Beta 1.2')
% figure
% imagesc(squeeze(nL2(:,54,:)))
% imagesc(squeeze(nL2(:,54,:)))
% title('Beta 2')
% figure
% imagesc(squeeze(nL5(:,54,:)))
% title('Beta 5')
% figure
% imagesc(squeeze(Axial(:,54,:)))
% colormap(gray)