function [ctxt,disto] = analyze_cb(cb)
%
% INPUT : 
% - CB : structure containing all information about the code-block to be
% processed. In particular, fields NBPS, COEFF and SIGN contain
% respectively the number of significant bit-planes, the magnitude of the
% wavelet coefficients and the sign of the wavelet coefficients.
%
% OUTPUTS : 
% - CTXT : matrix of size NBP x NCTXT x I where NBP is the number of
% bitplanes of the processed code-block, NCTXT is the number of contexts
% and I=2.
% It stores for each context and each bit-plane, the total number of bits
% (i=1) and the number of bits equal to '1' (i=2).
% - DISTO : vector of size NBP where NBP is the number of significant
% bit-planes of the processed code-block. Each element corresponds to the
% square error that would be obtained if the code-block was truncated at
% this bit-plane.

global nctxt;

ctxt = zeros(cb.nbps,nctxt,2);% 1 = nb total de bits, 2 = nb total de '1'
disto = zeros(cb.nbps,1);
sigstate = zeros(size(cb.coeff));
refstate = zeros(size(cb.coeff));

c_orig = cb.sign .* double(cb.coeff);

for bpno=1:cb.nbps % bpno 1 = msb
    
    nbp_discard = cb.nbps-bpno+1;
    
    %disto computation
    
    mask_AND = bitshift(uint16(65535),nbp_discard); % 65535 = maximum dynamic with UINT16
    mask_OR = bitshift(uint16(1),nbp_discard-1);
    if nbp_discard==1
        mask_OR=0;
    end
    m_trunc = bitor(bitand(cb.coeff,mask_AND),mask_OR);
    m_trunc(cb.coeff<=pow2(nbp_discard)-1)=0;
    
    c_trunc = cb.sign .* double(m_trunc);
    disto(bpno)=sum(sum((c_orig-c_trunc).^2,1),2);
    
    % Context computation
    
    for i=1:size(cb.coeff,2)
        for j=1:size(cb.coeff,1)
            b = bitand(cb.coeff(j,i),bitshift(uint16(1),nbp_discard-1));
            c = get_context(i,j,sigstate,refstate,cb.sbtype);
            if sigstate(j,i)
                refstate(j,i)=1;
            elseif b
                sigstate(j,i)=1;
            end
            ctxt(bpno,c,1)=ctxt(bpno,c,1)+1;
            if b
                ctxt(bpno,c,2)=ctxt(bpno,c,2)+1;
            end
        end
    end

end

return;