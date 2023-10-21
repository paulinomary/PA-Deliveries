function [rc rnc] = get_rate(ctxt,coeff,bpno,CDT)
% INPUTS
% - CTXT : matrix of size NBP x NCTXT x I where NBP is the number of
% bitplanes of the processed code-block, NCTXT is the number of contexts
% and I=2.
% It stores for each context and each bit-plane, the total number of bits
% (i=1) and the number of bits equal to '1' (i=2).
% - COEFF : coefficients of the processed code-block.
% - BPNO : bit-plane number
% - CDT : Context Distribution Table to be used for computing the rate.
% OUTPUTS
% - RC : Rate obtained when encoding the bit-plane based on the context
% distribution.
% - RNC : Rate obtained when encoding the bit-plane with a global
% distribution (no context). The global probability of getting a '1' is
% the last element of the given CDT. 

nbps = size(ctxt,1);
ctxt = permute(ctxt(bpno,:,:),[2 3 1]);
entropy = sum(-(log2(CDT(1:end-1))).*ctxt(:,2)-(log2((1-CDT(1:end-1)))).*(ctxt(:,1)-ctxt(:,2)));
nsign = numel(coeff(coeff<pow2(nbps-bpno+1) & coeff>=pow2(nbps-bpno)));
rc = entropy+nsign;
rnc = -(log2(CDT(end))).*sum(ctxt(:,2),1)-(log2(1-CDT(end))).*(sum(ctxt(:,1),1)-sum(ctxt(:,2),1))+nsign;

return;