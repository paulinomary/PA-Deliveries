function varargout = ndgrid(varargin)

% ndgrid for cells

if nargin == 0, 
    error('Not enough input arguments.'); 
end
if nargin == 1, 
    varargin = repmat(varargin, [1 max(nargout, 2)]);
end

nin = length(varargin);
nout = max(nargout,nin);

for i=length(varargin):-1:1,
%  varargin{i} = full(varargin{i}); % Make sure everything is full
  siz(i) = prod(size(varargin{i}));
end
if length(siz)<nout, siz = [siz ones(1,nout-length(siz))]; end

varargout = cell(1,nout);
for i=1:nout,
  x = varargin{i}(:); % Extract and reshape as a vector.
  s = siz; s(i) = []; % Remove i-th dimension
  x = reshape(x(:,ones(1,prod(s))),[length(x) s]); % Expand x
  varargout{i} = permute(x,[2:i 1 i+1:nout]); % Permute to i'th dimension
end
