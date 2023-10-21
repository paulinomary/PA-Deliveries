function si = imagesize(imagefile, varargin)

%imagesize   Size of image in file.
%
%   si = imagesize(imagefile[, format]) returns the size of the image stored in 
%   file 'imagefile'.  Format is the format of the file (see imfinfo), which may
%   be omited, in which case the format will be guessed, if possible.  The size
%   will always be returned has an integer row vector with length 3 (i.e., the 
%   number of color planes is always present, even if there is only one).
%
%   See imfinfo.
%
%Copyright (C) 2002 - Manuel Menezes de Sequeira (ISCTE).

info = imfinfo(imagefile, varargin{:});

if strcmp(info.ColorType, 'truecolor')
    si(3) = 3;
elseif strcmp(info.ColorType, 'grayscale')
    si(3) = 1;
elseif strcmp(info.ColorType, 'indexed')
    if isfield(info, 'Colormap')
        if mapisgray(info.Colormap)
            si(3) = 1;
        else
            si(3) = 3;
        end
    elseif isfield(info, 'ColorTable')
        if mapisgray(info.ColorTable)
            si(3) = 1;
        else
            si(3) = 3;
        end
    else
        error('indexed image: colormap not found');
    end
else
    error('unknown color type');
end
        
si(1) = info.Height;
si(2) = info.Width;