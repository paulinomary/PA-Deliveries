% This is a companion file to the book "Applied Signal Processing",
% by T.Dutoit and F. Marques, Springer 2008.
%
% It is supposed to be run cell-by-cell, using the cell mode of
% MATLAB 6 and later versions. Search for "what are cells" in the
% product help of MATLAB to know how to use them.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function calculates U(x) given a label l and a location
% in the labelmap L x,y,z
%
% considering a 10 neighbourhood

function [UX]=ASP_brain_segmentation_calcux_5class(L,Sx,Sy,Sz,label) 

% don't calculate UX at the boundaries where there are not 
% enough neighbours!

s=size(L);
NL=s(1);
NC=s(2);
NS=s(3);

%creation of a label matrix for every clique
disp('Computing UX: creating neighbourhood...');
NRight=[L(:,2:NC,:) zeros(NL,1,NS)];
NLeft=[zeros(NL,1,NS) L(:,1:(NC-1),:)];
NTop=[zeros(1,NC,NS); L(1:NL-1,:,:)];
NBottom=[L(2:NL,:,:); zeros(1,NC,NS)];
NFront=zeros(NL,NC,NS);
NFront(:,:,1:NS-1)=L(:,:,2:NS);
NBack=zeros(NL,NC,NS);
NBack(:,:,2:NS)=L(:,:,1:NS-1);
%Diagonal neihbourhood
NDTopLeft=[zeros(1,NC,NS); NLeft(1:NL-1,:,:)];
NDTopRight=[zeros(1,NC,NS); NRight(1:NL-1,:,:)];
NDBottomLeft=[NLeft(2:NL,:,:); zeros(1,NC,NS)];
NDBottomRight=[NRight(2:NL,:,:); zeros(1,NC,NS)];

%Following new cliques as presented by Shattuck
%atencio estic considerant el bck.... i canviare el seu valor...
TMP=zeros(NL,NC,NS);
TMP(find(NRight==label))=2;
TMP(find(NRight==label+1))=1;
TMP(find(NRight==label-1))=1;
TMP(find(NRight>label+1))=-1;
TMP(find(NRight<label-1))=-1;
clear NRight
UX=-TMP./Sx;

clear TMP
TMP=zeros(NL,NC,NS);
TMP(find(NLeft==label))=2;
TMP(find(NLeft==label+1))=1;
TMP(find(NLeft==label-1))=1;
TMP(find(NLeft>label+1))=-1;
TMP(find(NLeft<label-1))=-1;
clear NLeft
UX=UX-TMP./Sx;

clear TMP
TMP=zeros(NL,NC,NS);
TMP(find(NTop==label))=2;
TMP(find(NTop==label+1))=1;
TMP(find(NTop==label-1))=1;
TMP(find(NTop>label+1))=-1;
TMP(find(NTop<label-1))=-1;
clear NTop
UX=UX-TMP./Sy;

clear TMP
TMP=zeros(NL,NC,NS);
TMP(find(NBottom==label))=2;
TMP(find(NBottom==label+1))=1;
TMP(find(NBottom==label-1))=1;
TMP(find(NBottom>label+1))=-1;
TMP(find(NBottom<label-1))=-1;
clear NBottom
UX=UX-TMP./Sy;

clear TMP
TMP=zeros(NL,NC,NS);
TMP(find(NFront==label))=2;
TMP(find(NFront==label+1))=1;
TMP(find(NFront==label-1))=1;
TMP(find(NFront>label+1))=-1;
TMP(find(NFront<label-1))=-1;
clear NFront
UX=UX-TMP./Sz;

clear TMP
TMP=zeros(NL,NC,NS);
TMP(find(NBack==label))=2;
TMP(find(NBack==label+1))=1;
TMP(find(NBack==label-1))=1;
TMP(find(NBack>label+1))=-1;
TMP(find(NBack<label))=-1;
clear Nback
UX=UX-TMP./Sz;

clear TMP
TMP=zeros(NL,NC,NS);
TMP(find(NDTopLeft==label))=2;
TMP(find(NDTopLeft==label+1))=1;
TMP(find(NDTopLeft==label-1))=1;
TMP(find(NDTopLeft>label+1))=-1;
TMP(find(NDTopLeft<label-1))=-1;
clear NDTopLeft
UX=UX-TMP./sqrt(Sx^2+Sy^2);

clear TMP
TMP=zeros(NL,NC,NS);
TMP(find(NDTopRight==label))=2;
TMP(find(NDTopRight==label+1))=1;
TMP(find(NDTopRight==label-1))=1;
TMP(find(NDTopRight>label+1))=-1;
TMP(find(NDTopRight<label-1))=-1;
clear NDTopRight
UX=UX-TMP./sqrt(Sx^2+Sy^2);

clear TMP
TMP=zeros(NL,NC,NS);
TMP(find(NDBottomLeft==label))=2;
TMP(find(NDBottomLeft==label+1))=1;
TMP(find(NDBottomLeft==label-1))=1;
TMP(find(NDBottomLeft>label+1))=-1;
TMP(find(NDBottomLeft<label-1))=-1;
clear NDBottomLeft
UX=UX-TMP./sqrt(Sx^2+Sy^2);

clear TMP
TMP=zeros(NL,NC,NS);
TMP(find(NDBottomRight==label))=2;
TMP(find(NDBottomRight==label+1))=1;
TMP(find(NDBottomRight==label-1))=1;
TMP(find(NDBottomRight>label+1))=-1;
TMP(find(NDBottomRight<label-1))=-1;
clear NDBottomRight
UX=UX-TMP./sqrt(Sx^2+Sy^2);
clear TMP

%UX=UX;

UX(:,:,1)=ones(NL,NC);
UX(:,:,NS)=ones(NL,NC);
UX(1,:,:)=ones(NC,NS);
UX(NL,:,:)=ones(NC,NS);
UX(:,1,:)=ones(NL,NS);
UX(:,NC,:)=ones(NL,NS);

%save UX UX


