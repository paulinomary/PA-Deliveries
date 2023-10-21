function s = sequencesize(basename)

% sequencesize - Returns the size of the images images in a sequence.
%
%
%   s = sequencesize(basename)  Returns the size of the images in the
%   sequence stored in directory 'basename'.  Assumes all images in the directory 
%   have the same size.
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

d = dir(basename);
d = sort({d(3:end).name});

s = imagesize([basename '/' d{1}]);