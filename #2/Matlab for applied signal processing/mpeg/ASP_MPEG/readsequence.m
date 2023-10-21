function sequence = readsequence(basename, first, last, varargin)

% readsequence - Reads an image sequence.
%
%   Arguments:
%    basename - the name of the directory holding the sequence's images (all 
%               files in the directory must be images of the same type and size).
%    first - the first image of the sequence to read (defaults to 1).
%    last - the last image of the sequence to read (defaults to 1).
%    extention(s) - (optional) the extentions of the sequence images.
%           Note: first and last indexes refer only to image files
% 
%       Example: sequence = readsequence('.', 'jpg', 'bmp')
%
% Copyright (C) 2001 - Manuel Menezes de Sequeira (ISCTE).

if ~exist(basename, 'dir')
    error('basename must be an existing directory');
end
if(nargin < 2)
    first = 1;
end
if(nargin < 3)
    last = 1;
end

if nargin > 3
    d = dir(strcat(basename, '/*.', varargin{1}));
    for i = 2 : length(varargin)
        d = cat(1, d, dir(strcat(basename, '/*.', varargin{i})));
    end
    d = sort({d.name});
else
    d = dir(basename);
    d = sort({d(3:end).name});
end

if first > length(d)
    sequence = zeros(1, 1, 1, 0);
    return;
end

first = max(1, first);
last = min(size(d, 2), last);

if first > last
    sequence = zeros(1, 1, 1, 0);
    return;
end

for i = 1 : min(size(d, 2), last) - max(1, first) + 1,
    sequence(:, :, :, i) = readimage([basename '/' d{i + max(1, first) - 1}]);
end