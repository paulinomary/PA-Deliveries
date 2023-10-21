function i = readimage(file)

% readimage - Reads an image from a file.
%
%
%   i = readimage(file)  Substitute for imread.  Does not use colormaps!
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

[i, map] = imread(file);

if(~isempty(map))
    if(mapisgray(map))
        if mapisbinary(map)
            i = logical(uint8(i));
            return
        end
        if(isa(i, 'uint8'))
            i = uint8(ind2gray(i, map)*255);
        else
            i = ind2gray(i, map);
        end
    else
        if(isa(i, 'uint8'))
            i = uint8(ind2rgb(i, map)*255);
        else    
            i = ind2rgb(i, map);
        end
    end
end

if ndims(i) < 3
    v = unique(i);
    if(length(v) == 2)
        % It is binary?
        if v(1) == 0 & (v(2) == 255  | v(2) == 1)
            i = i ~= 0;
        end
    end
end
            