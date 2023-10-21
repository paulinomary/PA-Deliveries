function D = get_disto(disto,bpno)
% INPUTS : 
% - DISTO : vector of size NBP where NBP is the number of significant
% bit-planes of the processed code-block. Each element corresponds to the
% square error that would be obtained if the code-block was truncated at
% this bit-plane.
% - BPNO : the bit-plane for which we want to compute the distortion
% reduction brought its decoding.
%
% OUTPUT : D, the distortion reduction brought by the decoding of a given
% bit-plane.

if bpno<numel(disto)
    D = disto(bpno)-disto(bpno+1);
else
    D = disto(end);
end

return;