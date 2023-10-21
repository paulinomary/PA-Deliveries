% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.

close all
clear all
warning off MATLAB:divideByZero

%%
% Gui to open a the 3D brain data file
disp('')
disp('Please select a 3D brain image file')
[bfilename, bpathname] = uigetfile('*.*', 'Open a .mat file');
filenameMRI=strcat(bpathname,bfilename);
load(filenameMRI);
%%
%Let's first visualize the 3D image data
%Image at time 1
   load Image1_voi.mat
   figure(1); clf;
   subplot(1,3,1); imagesc(squeeze(Axial(90,:,:)));  
     iptsetpref('ImshowBorder', 'tight');
    colormap(gray);
     title('Coronal view'); hold on; 
   subplot(1,3,2); imagesc(squeeze(Axial(:,65,:)));  
    iptsetpref('ImshowBorder', 'tight'); 
    title('Axial view'); hold on;
   subplot(1,3,3); imagesc((squeeze(Axial(:,:,60)))'); 
     axis fill;
     iptsetpref('ImshowBorder', 'tight');
     colorbar; title('Sagittal view');
%Image at time 2
   load Image2_voi.mat
   figure(2); clf;
   subplot(2,3,1); imagesc(squeeze(Axial(90,:,:)));  
     iptsetpref('ImshowBorder', 'tight');
     colormap(gray);
     title('Coronal view'); hold on; 
   subplot(2,3,2); imagesc(squeeze(Axial(:,65,:)));  
    iptsetpref('ImshowBorder', 'tight'); 
    title('Axial view'); hold on;
   subplot(2,3,3); imagesc((squeeze(Axial(:,:,60)))'); 
     axis fill;
     iptsetpref('ImshowBorder', 'tight');
     colorbar; title('Sagittal view');
% Image volume
   load ImageVolume.mat
   figure(3); clf;
   phandles = contourslice(Axial,[],[45,64,78,90,100,110,120,127,135,142,150,160,175,190,200,210,220,228],[60,65,70,75,80,85],7);
   view(-169,-54); axis tight; set(phandles,'LineWidth',2); colormap(gray); title('3D MRI');

%%
% Let's have a look on the image histogram
load Image1_voi.mat
loffset=1; % Avoid the pic of the histogram at 0 gray-level due to background
roffset=0; 
[H,X,NbrBin,BinSize]=ASP_brain_segmentation_histogram(Axial,loffset,roffset);
figure(4);plot(X,H);
title('Intensity Image 1 Histogram')
xlabel('Gray values')
ylabel('Number of elements')
axis([1 600 0 8000]) %Axial1

load Image2_voi.mat
loffset=1; % Avoid the pic of the histogram at 0 gray-level due to background
roffset=0; 
[H,X,NbrBin,BinSize]=ASP_brain_segmentation_histogram(Axial,loffset,roffset);
figure(5);plot(X,H);
title('Intensity Image 2 Histogram')
xlabel('Gray values')
ylabel('Number of elements')
axis([1 150 0 40000]) %Axial2

%%
% To fit just a part of the histogram set manually loffset (left offset) and/or
% roffset (right offset)
% gui to open a the 3D brain data file
disp('')
disp('Please select a 3D brain image file')
[bfilename, bpathname] = uigetfile('*.*', 'Open a .mat file');
filenameMRI=strcat(bpathname,bfilename);
load(filenameMRI);

loffset=1; % Avoid the pic of the histogram at 0 gray-level due to background
roffset=0; 
[H,X,NbrBin,BinSize]=ASP_brain_segmentation_histogram(Axial,loffset,roffset);
figure(3);plot(X,H);
title('Intensity Image Histogram')
xlabel('Gray values')
ylabel('Number of elements')

% get loffset and roffset
disp(' ')
disp('Click the UPPER and LOWER boundary of the histogram')
disp('(First click on the LEFT boundary with the LEFT mouse button,')
disp('then click on the RIGHT boundary with the RIGHT mouse button)');
disp(' ');
disp('press any key when ready...');
pause
lims=getpts;
disp(' ')
loffset=round(lims(1))
roffset=round(NbrBin-lims(2))
[H,X,NbrBin,BinSize]=ASP_brain_segmentation_histogram(Axial,loffset,roffset);
figure(3);plot(X,H);
title('Intensity Image Histogram')
xlabel('Gray values')
ylabel('Number of elements')

%%
% Get initial values for the EM algorithm
% Get initial mu's by clicking on the image histogram
% disp(' ')
% disp('Please click on the MEAN values of the Gaussians to fit.')
% disp('Only the x-axis values will be used.')
% disp('If you make an error remove the point with backspace.')
% disp('Click THE LAST MEAN VALUE with the RIGHT MOUSE BUTTON')
% disp(' ')
% disp('press any key when ready');
% pause
% disp('please wait...')
loffset=1; % Avoid the pic of the histogram at 0 gray-level due to background
roffset=0; 
[H,X,NbrBin,BinSize]=ASP_brain_segmentation_histogram(Axial,loffset,roffset);
NbrSamples=sum(sum(H));
nH=H*1/NbrSamples;
% figure(3);
% plot(X,H);
% title('click mu-values')
% xlabel('Gray values')
% mus=sortrows(getpts); % To pick manually the mean values

loffset=1; % Avoid the pic of the histogram at 0 gray-level due to background
roffset=0; 
%mus = [130 220 330]'; % image 1 but wrong convergence
mus = [115 345 410]'; % image 1
% mus = [43 100 123]'; % image 2
s=size(mus);
NbrComp=s(1);

for i=1:NbrComp
	sig(i,1)=5.0; % Small value different from zero
end

for i=1:NbrComp
	%pw(i,1)=input('probablity=');
    pw(i,1)=1/NbrComp; % Initially, the assumption that all classes are equiprobable is done
end

% compose Theta
Theta=cat(2,mus,sig,pw)
% ---
 
% starting calculations
disp(' ');
disp('Now starting EM algorithm...');
 
[S,e,Theta,fitper,G]=ASP_brain_segmentation_segmentem(Axial, loffset, roffset, Theta, BinSize);
 
disp('Estimated means and variance:');
disp('Mu');
ul=Theta(:,1)
disp('Sigma');
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

%%
% Reading GMM-EM classification
% Open label map 
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

% how many iterations?
% nit=input('How many iterations must be done for the HMM: ');
nit=5;
beta=0.6;
nLold=nL;
Change=zeros(nit,1);
for ii=1:nit
	disp('Now doing the MAP part');
	nL=ASP_brain_segmentation_daloop(nL,Axial,ul,sl,NbrComp,beta);
    disp('New map nL calculated');
    disp('Computing number of voxels that changed');
    Diff=(nLold-nL);
    Change(ii)=length(find(Diff))./length(find(nL))*100
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
%%
%Compute probability maps
load Image1_voi.mat
%load Image2_voi.mat

s=size(Axial);

% Open HMRF-EM classification 
fid=fopen('Image1_voi_nL.raw','r','l');
nL=fread(fid,s(1)*s(2)*s(3),'ushort');
nL=reshape(nL,s(1),s(2),s(3));

%filenameMRI,ul and sl are in memory from other cells
ASP_brain_segmentation_ProbMaps(Axial, nL, ul, sl, NbrComp,filenameMRI(1:end-4));

% Smooth probability maps = Concentration maps
Ccsf_1=ASP_brain_segmentation_smooth_3DGauss('Image1_voi_Pcsf.raw',s(1),s(2),s(3),'float',7,11);
Cgm_1=ASP_brain_segmentation_smooth_3DGauss('Image1_voi_Pgm.raw',s(1),s(2),s(3),'float',7,11);
Ccsf_2=ASP_brain_segmentation_smooth_3DGauss('Image2_voi_Pcsf.raw',s(1),s(2),s(3),'float',7,11);
Cgm_2=ASP_brain_segmentation_smooth_3DGauss('Image2_voi_Pgm.raw',s(1),s(2),s(3),'float',7,11);

% Visualize
figure;
   subplot(1,3,1); imagesc(squeeze(Cgm_1(90,:,:)));  
     iptsetpref('ImshowBorder', 'tight');
      title('Coronal view'); hold on; 
   subplot(1,3,2); imagesc(squeeze(Cgm_1(:,65,:)));  
    iptsetpref('ImshowBorder', 'tight'); 
    title('Axial view'); hold on;
   subplot(1,3,3); imagesc((squeeze(Cgm_1(:,:,60)))'); 
     axis fill;
     iptsetpref('ImshowBorder', 'tight');
     colorbar; title('Sagittal view');
     
 figure;
   subplot(1,3,1); imagesc(squeeze(Ccsf_1(90,:,:)));  
     iptsetpref('ImshowBorder', 'tight');
      title('Coronal view'); hold on; 
   subplot(1,3,2); imagesc(squeeze(Ccsf_1(:,65,:)));  
    iptsetpref('ImshowBorder', 'tight'); 
    title('Axial view'); hold on;
   subplot(1,3,3); imagesc((squeeze(Ccsf_1(:,:,60)))'); 
     axis fill;
     iptsetpref('ImshowBorder', 'tight');
     colorbar; title('Sagittal view');

%%
% Computing the regions where a significant difference between time 1 and
% time 2 exists
D12=abs(Ccsf_1-Ccsf_2).*abs(Cgm_1-Cgm_2);
figure;imagesc(squeeze(D12(:,99,:))); colorbar;
axis([2 121 40 186]);title('Tissue degeneration between time 1 and 2');

figure;imagesc(squeeze(D12(84,:,:)));colorbar;
title('Tissue degeneration between time 1 and 2');

% Region of interest: threshold of the degeneration map (at 500 for
% instance)
volume_roi_1=sum(Cgm_1(D12>500))
volume_roi_2=sum(Cgm_2(D12>500))
degeneration_roi=(volume_roi_1-volume_roi_2)/volume_roi_1*100

% Region of interest (entorhinal cortex)
ROD = D12(45:115,70:120,70:121);
% Percentage of GM in this region
volume_roi_1=sum(Cgm_1(ROD~=0))
volume_roi_2=sum(Cgm_2(ROD~=0))
degeneration_roi=(volume_roi_1-volume_roi_2)/volume_roi_1*100

%%
warning on MATLAB:divideByZero
