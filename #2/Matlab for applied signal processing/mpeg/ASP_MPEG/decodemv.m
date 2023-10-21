function varargout = decodemv(codedmessage, imsize, bs, varargin)

% DECODEMV Motion vectors encoder.
%
%   [MVV, MHV, CODEDMESSAGE] = DECODEMV(CODEDMESSAGE, IMSIZE, BS) decodes
%   message and returns the motion vertical and horizontal vectores values
%   (MVV and MHV) encoded in CODEDMESSAGE. IMSIZE is the size of the
%   decoded image and must be multiple of BS. BS is the size of one block,
%   typically [16 16] or [8 8].
%   
%   [MVV, MHV, CODEDMESSAGE] = DECODEMV(CODEDMESSAGE, IMSIZE, BS, TYPE)
%   does the same above, but with type specified. TYPE must be 'motion'.
%
%   [MVV1, MHV1, MVV2, MHV2, CODEDMESSAGE] = ...
%      DECODEMV(CODEDMESSAGE, IMSIZE, BS, TYPE, MASK1, MASK2) in this case
%   CODEDMESSAGE contains four motion vectores, two vertical and two
%   horizontal. TYPE must be 'motion2'. MASK1 and MASK2 are the same used
%   at the encoder and determin the position of the vector values.
%   
% See ENCODEMV

% Calculates the motion vertical and horizontal vector values length.
blksize = imsize ./ bs;
mvlen = prod(blksize);

if nargin > 3
    type = varargin{1};
else
    type = 'motion';
end

switch type
    case 'motion'
        load('./decoders/mv/mvvdecoder.mat');
        [mvvindex codedmessage] = decodemessage(codedmessage, mvvdecoder, mvlen);
        clear mvvdecoder;
        load('./decoders/mv/mhvdecoder.mat');
        [mhvindex codedmessage] = decodemessage(codedmessage, mhvdecoder, mvlen);
        clear mhvdecoder;
        varargout{1} = reshape(mvvindex', blksize) - 16;
        varargout{2} = reshape(mhvindex', blksize) - 16;
        varargout{3} = codedmessage;
        
    case 'motion2'
        mask1 = varargin{2}; 
        mask2 = varargin{3};
        
        nBits = ceil(log2(prod(blksize)));
        mvlen1 = bin2dec(codedmessage(1, 1 : nBits));
        mvlen2 = bin2dec(codedmessage(1, nBits + 1 : 2 * nBits));
        
        load('./decoders/mv/mvvdecoder.mat');
        [mvvindex1 codedmessage] = decodemessage(codedmessage(1, 2 * nBits + 1 : end), mvvdecoder, mvlen1);
        load('./decoders/mv/mhvdecoder.mat');
        [mhvindex1 codedmessage] = decodemessage(codedmessage, mhvdecoder, mvlen1);
        [mvvindex2 codedmessage] = decodemessage(codedmessage, mvvdecoder, mvlen2);
        clear mvvdecoder;
        [mhvindex2 codedmessage] = decodemessage(codedmessage, mhvdecoder, mvlen2);        
        clear mhvdecoder;
        
        mvv1 = zeros(blksize);
        mvv1(mask1) = mvvindex1 - 16;
        mhv1 = zeros(blksize);
        mhv1(mask1) = mhvindex1 - 16;
        
        mvv2 = zeros(blksize);
        mvv2(mask2) = mvvindex2' - 16;
        mhv2 = zeros(blksize);
        mhv2(mask2) = mhvindex2' - 16;
        
        varargout{1} = mvv1;
        varargout{2} = mhv1;
        varargout{3} = mvv2;
        varargout{4} = mhv2;
        varargout{5} = codedmessage;
        
    otherwise
        error('Type not supported!! Check help decodemv!!');
end