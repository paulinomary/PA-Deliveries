function codedmessage = encodemv(mvv, mhv, varargin)

% ENCODEMV Motion vectors encoder.
%
%   CODEDMESSAGE = ENCODEMV(MVV, MHV) encodes vertical and horizontal
%   motion vectors values, and returns a char row vector with the coded
%   message.
%   
%   CODEDMESSAGE = ENCODEMV(MVV, MHV, TYPE) does the same as above, but
%   type is specified, and must be 'motion'.
%
%   CODEDMESSAGE = ENCODEMV(MVV1, MHV1, TYPE, MVV2, MVV2, MASK1, MASK2) in
%   this case it encodes two motion vectors, MVV1/MHV1, and MVV2/MHV2. TYPE
%   must be 'motion2'. MASK1 and MASK2 select witch components of the
%   motion vectors values are selected.
%   
% See DECODEMV

if nargin > 2
    type = varargin{1};
else
    type = 'motion';
end

switch type
    case 'motion'
        load('./codes/mvcodes.mat');
        codedmvv = encodemessage(mvv(:)' + 16, mvvcode);
        codedmhv = encodemessage(mhv(:)' + 16, mhvcode);
        codedmessage = strcat(codedmvv, codedmhv);
  
    case 'motion2'
        if nargin < 7
            error('Not enough input arguments!!');
        else
            load('./codes/mvcodes.mat');
            mvv2 = varargin{2};
            mhv2 = varargin{3};
            mask1 = varargin{4}; 
            mask2 = varargin{5};
            
            % Calculate number of bits used to encode motion vectores
            % length.
            nbits = ceil(log2(prod(size(mvv))));
            
            % Select motion vectores used and discard thouse that aren't
            % used.
            mvv = mvv(mask1);
            mhv = mhv(mask1);
            mvv2 = mvv2(mask2);
            mhv2 = mhv2(mask2);
            
            % Calculate motion vectors length.
            mvlen1 = length(mvv);
            mvlen2 = length(mvv2);
            
            codedmessage = strcat(dec2bin(mvlen1, nbits), ...
                dec2bin(mvlen2, nbits));
            
            if mvlen1 > 0
                mvvcm1 = encodemessage(mvv' + 16, mvvcode);
                mhvcm1 = encodemessage(mhv' + 16, mhvcode);
                codedmessage = strcat(codedmessage, mvvcm1, mhvcm1);
            end
            
            if mvlen2 > 0
                mvvcm2 = encodemessage(mvv2' + 16, mvvcode);
                mhvcm2 = encodemessage(mhv2' + 16, mhvcode);
                codedmessage = strcat(codedmessage, mvvcm2, mhvcm2);
            end
        end
        
    otherwise
        error('Type not supported!! Check help encodemv!!');
end