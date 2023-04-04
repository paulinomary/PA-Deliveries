

function Callbacks_uniform_quantization_GUI25(f,C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%SENSE COMPUTER AND SET FILE DELIMITER
switch(computer)				
    case 'MACI64',	char= '/';
    case 'GLNX86',  char='/';
    case 'PCWIN',	char= '\';
    case 'PCWIN64', char='\';
    case 'GLNXA64', char='/';
end

x=C{1,1};
y=C{1,2};
a=C{1,3};
b=C{1,4};
u=C{1,5};
v=C{1,6};
m=C{1,7};
n=C{1,8};
lengthbutton=C{1,9};
widthbutton=C{1,10};
enterType=C{1,11};
enterString=C{1,12};
enterLabel=C{1,13};
noPanels=C{1,14};
noGraphicPanels=C{1,15};
noButtons=C{1,16};
labelDist=C{1,17};%distance that the label is below the button
noTitles=C{1,18};
buttonTextSize=C{1,19};
labelTextSize=C{1,20};
textboxFont=C{1,21};
textboxString=C{1,22};
textboxWeight=C{1,23};
textboxAngle=C{1,24};
labelHeight=C{1,25};
fileName=C{1,26};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %PANELS
% for j=0:noPanels-1
% uipanel('Parent',f,...
% 'Units','Normalized',...
% 'Position',[x(1+4*j) y(1+4*j) x(2+4*j)-x(1+4*j) y(3+4*j)-y(2+4*j)]);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GRAPHIC PANELS
for i=0:noGraphicPanels-1
switch (i+1)
case 1
graphicPanel1 = axes('parent',f,...
'Units','Normalized',...
'Position',[a(1+4*i) b(1+4*i) a(2+4*i)-a(1+4*i) b(3+4*i)-b(2+4*i)],...
'GridLineStyle','--');
case 2
graphicPanel2 = axes('parent',f,...
'Units','Normalized',...
'Position',[a(1+4*i) b(1+4*i) a(2+4*i)-a(1+4*i) b(3+4*i)-b(2+4*i)],...
'GridLineStyle','--');
case 3
graphicPanel3 = axes('parent',f,...
'Units','Normalized',...
'Position',[a(1+4*i) b(1+4*i) a(2+4*i)-a(1+4*i) b(3+4*i)-b(2+4*i)],...
'GridLineStyle','--');
case 4
graphicPanel4 = axes('parent',f,...
'Units','Normalized',...
'Position',[a(1+4*i) b(1+4*i) a(2+4*i)-a(1+4*i) b(3+4*i)-b(2+4*i)],...
'GridLineStyle','--');
case 5
graphicPanel5 = axes('parent',f,...
'Units','Normalized',...
'Position',[a(1+4*i) b(1+4*i) a(2+4*i)-a(1+4*i) b(3+4*i)-b(2+4*i)],...
'GridLineStyle','--');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TITLE BOXES
for k=0:noTitles-1
switch (k+1)
case 1
titleBox1 = uicontrol('parent',f,...
'Units','Normalized',...
'Position',[u(1+4*k) v(1+4*k) u(2+4*k)-u(1+4*k) v(3+4*k)-v(2+4*k)],...
'Style','text',...
'FontSize',textboxFont{k+1},...
'String',textboxString(k+1),...
'FontWeight',textboxWeight{k+1},...
'FontAngle',textboxAngle{k+1});
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUTTONS
for i=0:(noButtons-1)
enterColor='w';
if strcmp(enterType{i+1},'pushbutton')==1 ||strcmp(enterType{i+1},'text')==1
enterColor='default';
end
if (strcmp(enterLabel{1,(i+1)},'')==0 &&...
        strcmp(enterLabel{1,(i+1)},'...')==0) %i.e. there is a label
%creating a label for some buttons
uicontrol('Parent',f,...
'Units','Normalized',...
'Position',[m(1+2*i) n(1+2*i)-labelDist-labelHeight(i+1) ...
(m(2+2*i)-m(1+2*i)) labelHeight(i+1)],...
'Style','text',...
'String',enterLabel{i+1},...
'FontSize', labelTextSize(i+1),...
'HorizontalAlignment','center');
end
switch (i+1)
case 1
button1=uicontrol('Parent',f,...
'Units','Normalized',...
'Position',[m(1+2*i) n(1+2*i) (m(2+2*i)-m(1+2*i)) (n(2+2*i)-n(1+2*i))],...
'Style',enterType{i+1},...
'String',enterString{i+1},...
'FontSize', buttonTextSize(1+i),...
'BackgroundColor',enterColor,...
'HorizontalAlignment','center',...
'Callback',@button1Callback);
case 2
button2=uicontrol('Parent',f,...
'Units','Normalized',...
'Position',[m(1+2*i) n(1+2*i) (m(2+2*i)-m(1+2*i)) (n(2+2*i)-n(1+2*i))],...
'Style',enterType{i+1},...
'String',enterString{i+1},...
'FontSize', buttonTextSize(1+i),...
'BackgroundColor',enterColor,...
'HorizontalAlignment','center',...
'Callback',@button2Callback);
case 3
button3=uicontrol('Parent',f,...
'Units','Normalized',...
'Position',[m(1+2*i) n(1+2*i) (m(2+2*i)-m(1+2*i)) (n(2+2*i)-n(1+2*i))],...
'Style',enterType{i+1},...
'String',enterString{i+1},...
'FontSize', buttonTextSize(1+i),...
'BackgroundColor',enterColor,...
'HorizontalAlignment','center',...
'Callback',@button3Callback);
case 4
button4=uicontrol('Parent',f,...
'Units','Normalized',...
'Position',[m(1+2*i) n(1+2*i) (m(2+2*i)-m(1+2*i)) (n(2+2*i)-n(1+2*i))],...
'Style',enterType{i+1},...
'String',enterString{i+1},...
'FontSize', buttonTextSize(1+i),...
'BackgroundColor',enterColor,...
'HorizontalAlignment','center',...
'Callback',@button4Callback);
case 5
button5=uicontrol('Parent',f,...
'Units','Normalized',...
'Position',[m(1+2*i) n(1+2*i) (m(2+2*i)-m(1+2*i)) (n(2+2*i)-n(1+2*i))],...
'Style',enterType{i+1},...
'String',enterString{i+1},...
'FontSize', buttonTextSize(1+i),...
'BackgroundColor',enterColor,...
'HorizontalAlignment','center',...
'Callback',@button5Callback);
case 6
button6=uicontrol('Parent',f,...
'Units','Normalized',...
'Position',[m(1+2*i) n(1+2*i) (m(2+2*i)-m(1+2*i)) (n(2+2*i)-n(1+2*i))],...
'Style',enterType{i+1},...
'String',enterString{i+1},...
'FontSize', buttonTextSize(1+i),...
'BackgroundColor',enterColor,...
'HorizontalAlignment','center',...
'Callback',@button6Callback);
case 7
button7=uicontrol('Parent',f,...
'Units','Normalized',...
'Position',[m(1+2*i) n(1+2*i) (m(2+2*i)-m(1+2*i)) (n(2+2*i)-n(1+2*i))],...
'Style',enterType{i+1},...
'String',enterString{i+1},...
'FontSize', buttonTextSize(1+i),...
'BackgroundColor',enterColor,...
'HorizontalAlignment','center',...
'Callback',@button7Callback);
case 8
button8=uicontrol('Parent',f,...
'Units','Normalized',...
'Position',[m(1+2*i) n(1+2*i) (m(2+2*i)-m(1+2*i)) (n(2+2*i)-n(1+2*i))],...
'Style',enterType{i+1},...
'String',enterString{i+1},...
'FontSize', buttonTextSize(1+i),...
'BackgroundColor',enterColor,...
'HorizontalAlignment','center',...
'Callback',@button8Callback);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER CODE FOR THE VARIABLES AND CALLBACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Variables
    curr_file=1;
    fs=8000;
    directory_name='abcd';
    wav_file_names='abce';
    fin_path='filename';
    fname='output';
    nsamp=1;
    ss=1300;
    es=18800;
    nbits=10;
    nbins=51;
    xin=[];
    xinn=[];
    xinq=[];
    einq=[];
    Nfft=1024;
    Nwin=512;
    filein='input';

% Name the GUI
    set(f,'Name','uniform_quantize');

% CALLBACKS
% Callback for button1 -- Get Speech Files Directory
 function button1Callback(h,eventdata)
%% ***** modified below **************************************************************************
     if isempty(getpref('SpeechApps'))
         url = sprintf('%s%s',...
             'http://www.mathworks.com/matlabcentral/fileexchange/',...
             'submissions/45293/v/1/download/zip');
         [saveloc, ~, ~] = fileparts(pwd); %save to one level up from current folder
         % Create a waitbar during download
         h = waitbar(0.35,'This may take several minutes...',...
             'Name','Downloading Speech Files...');
         % Download the zipped file
         [filestr,status] = urlwrite(url,[saveloc filesep 'speech_files.zip'],...
             'Timeout',10);
         if status
             delete(h);
             hh1= helpdlg('Downloaded. Select a location to UNZIP the speech files.');
             uiwait(hh1);
             unziploc = uigetdir(saveloc,'Select a location to unzip the speech files');
             h2 = waitbar(0.2,'This may take a minute...',...
                 'Name','Unzipping the Speech Files to Location Selected...');
             unzip(filestr,unziploc);
             delete(h2)
             addpref('SpeechApps','path',unziploc);
             hh2= helpdlg('Ready. Select the speech_files folder in the next window');
             uiwait(hh2);
         else
             warndlg('No Internet Connection to MATLAB Central!');
         end
         
     else
     end
     directory_name = uigetdir(getpref('SpeechApps','path'));
%% ***** modified above *******************************************
     A=strvcat(strcat((directory_name),char,'*.wav'));
     struct_filenames=dir(A);
     wav_file_names={struct_filenames.name};
     set(button2,'String',wav_file_names);
     set(button2,'val',1);
     
% once the popupmenu/drop down menu is created, by default, the first
% selection from the popupmenu/drop down menu id not called
    indexOfDrpDwnMenu=1;
    
% by default first option from the popupmenu/dropdown menu will be loaded
    [curr_file,fs]=loadSelection(directory_name,wav_file_names,indexOfDrpDwnMenu);
 end

% Callback for button2 -- Choose speech file for play and plot
 function button2Callback(h,eventdata)
     indexOfDrpDwnMenu=get(button2,'val');
     [curr_file,fs]=loadSelection(directory_name,wav_file_names,indexOfDrpDwnMenu);
 end

%*************************************************************************
% function -- load selection from designated directory and file
%
function [curr_file,fs]=loadSelection(directory_name,wav_file_names,...
    indexOfDrpDwnMenu);
%
% read in speech/audio file
% fin_path is the complete path of the .wav file that is selected
    fin_path=strcat(directory_name,'/',strvcat(wav_file_names(indexOfDrpDwnMenu)));
    
% clear speech/audio file
    clear curr_file;
    
% read in speech/audio signal into curr_file; sampling rate is fs 
    [curr_file,fs]=audioread(fin_path);
    
% scale input to range -1 to +1
    xmax=max(abs(curr_file));
    xin=curr_file/xmax;
    
% create title information with file, sampling rate, number of samples
    fname=wav_file_names(indexOfDrpDwnMenu);
    FS=num2str(fs);
    nsamp=num2str(length(curr_file));
    file_info_string=strcat('  file: ',fname,', fs: ',FS,' Hz, nsamp:',nsamp);
    
% load filename (fname) from cell array
    fname=wav_file_names{indexOfDrpDwnMenu};
end

% Callback for button3 -- nbits: number of bits for quantizing samples
 function button3Callback(h,eventdata)
     nbits=str2num(get(button3,'string'));
     if ~((nbits >= 2 && nbits <= 16))
        waitfor(errordlg('nbits must be a positive integer between 2 and 16'))
        return;
     end
     nbits=round(nbits);
     set(button3,'string',num2str(nbits));
     button4Callback(h,eventdata);
 end

% Callback for button4 -- Run Uniform Quantize
 function button4Callback(h,eventdata)
     
% check editable buttons for changes
    % button3Callback(h,eventdata);
    
    ss=1;
    es=length(xin);
    
% run setup_uniform_quantize
    [xinn,xinq,einq]=setup_uniform_quantize(xin,fs,ss,es,nbits,nbins,fname);
 end

%**********************************************************************
    function [xinn,xinq,einq]=setup_uniform_quantize(xin,fs,ss,es,nbits,nbins,filein)
%
% quantize a speech signal using a uniform quantizer with B-
% bits per sample using truncation and saturation

% save speech samples from ss to es in file xinn
    xinn=xin(ss:es);
    
% normalize input so that peak is in the range [-1 1]
    xinn=xinn/max(max(xinn),-min(xinn));
    
% create synthetic ramp to test code -- used when nbits=0
    if (nbits == 0)
        xinn=-1:2/(es-1):1;
        nbits=4;
        set(button3,'string',num2str(nbits));
    end
    
% clear graphics Panel 5
        reset(graphicPanel5);
        axes(graphicPanel5);
        cla;
        
% plot original speech signal in graphics Panel 5
    stitle=sprintf('filein:%s, ss,es:%d %d, unquantized, xmin,xmax:%7.4f %7.4f',...
        filein,ss,es,min(xinn),max(xinn));
    n=ss-1:es-1;
    xpp=['Time in Samples; fs=',num2str(fs),' samples/second'];
    plot(n,xinn,'k');xlabel(xpp);ylabel('Value');
    legend('original speech file');
    grid on, axis tight;
    
% obtain estimate of long-time average power spectrum
    Nwin=512;
    Nfft=1024;
    [P1,F,R,T]=pd_spect_U(xinn,fs,Nfft,Nwin);
    
% clear graphics Panel 2
    reset(graphicPanel2);
    axes(graphicPanel2);
    cla;
    
% display the original signal power density spectrum in graphics Panel 2
    plot(F,10*log10(P1),'r','LineWidth',2),grid on,ylabel('Log Magnitude in dB'),...
        xlabel('Frequency in Hz'); hold on;
    stitle=sprintf('filein:%s, original/4-10 bit quantization, Nwin:%d',...
        filein,Nwin);
    
% quantize using nbits per sample and plot log spectrum of the 
% unquantized speech along with the log spectrum of the quantization signal
        delta=2/2^nbits;
        xinq=fxquant(xinn,nbits,'trunc','sat')+delta/2;
        
% form error signal
        einq=xinn-xinq;
        
% clear graphics Panel 4
        reset(graphicPanel4);
        axes(graphicPanel4);
        cla;  
        
% plot quantized signal in graphics Panel 4
        stitle=sprintf('filein:%s, quantized signal (B: %d bit quantization), Nwin:%d',...
        filein,nbits,Nwin);
        plot(n,xinq);xlabel(xpp);ylabel('Value');legend('quantized signal');
        grid on;axis tight;

% clear graphics Panel 3
        reset(graphicPanel3);
        axes(graphicPanel3);
        cla;
        
% plot quantization error signal in graphics Panel 3
        stitle=sprintf('filein:%s, error signal (B: %d bit quantization), Nwin:%d',...
        filein,nbits,Nwin);
        plot(n,einq);xlabel(xpp);ylabel('Value');legend('quantization error signal');
        grid on;axis tight;
            
% clear graphics Panel 1
            reset(graphicPanel1);
            axes(graphicPanel1);
            cla;
            
% plot error signal histogram in graphics Panel 1
        stitle=sprintf('filein:%s, error historgram at nbits:%d, with nbins:%d',...
            filein,nbits,nbins);
            hist(einq,nbins),grid on;legend('error histogram');
            % axis tight;
            
% obtain estimate of long-time average power spectrum
        hold on;
        [P2,F,R,T]=pd_spect_U(einq,fs,Nfft,Nwin);
        
% choose graphics Panel 2 for next plot
        axes(graphicPanel2); 
        
% plot error signal spectrum in graphics Panel 2
        plot(F,10*log10(P2),'b','LineWidth',2),grid on;
        G(1:length(F))=-6*nbits-4.77;
        
% plot theoretical error power spectrum as dashed line
        plot(F,G,'k--','LineWidth',2);
        legend('signal power spectrum','error power spectrum','theoretical error power spectrum');
        
%  compute signal-to-noise ratio (snr)
        [s_n_r,e]=snr(xinq,xinn);
        fprintf('uniform quantization, nbits:%d, snr:%7.2f \n',nbits,s_n_r);
        
% display fname and signal processing parameters in titleBox1
    stitle=sprintf(' file: %s, ss/es: %d  %d, nbits: %d',filein,ss,es,nbits); 
    stitle1=strcat('Uniform Quantization -- ',stitle);
    set(titleBox1,'string',stitle1);
    set(titleBox1,'FontSize',20);     
end

% callback for button 5 -- play original speech signal
function button5Callback(h,eventdata)
	soundsc(xinn,fs);
end
 
% callback for button 6 -- play quantized speech signal
function button6Callback(h,eventdata)
	soundsc(xinq,fs);
end
 
% callback for button 7 -- play error signal
function button7Callback(h,eventdata)
	soundsc(einq,fs);
 end

% Callback for button8 -- close GUI
 function button8Callback(h,eventdata)
     fclose('all');
     close(gcf);
 end
end