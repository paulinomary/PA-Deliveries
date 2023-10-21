function [type, index] = shufflestream(pP, pI, nframes)

type(1:nframes) = 'B';
type(1:pP:end) = 'P';
type(1:pP*pI:end) = 'I';

index = [1 : nframes] - 1;
index(1) = 1;

if type(end) == 'B'
    nI = ceil(nframes / (pP * pI));
    nP = ceil(nframes / pP) - nI;
    nB = (nP + nI - 1) * (pP - 1);
    nBleft = nframes - nB - nP - nI;
    index(2 : pP : end - nBleft) = index(2 : pP : end - nBleft) + pP;
    index(end - nBleft + 1) = nframes;

    if (nP - nI + 1) == pI - 1
        type(end) = 'I';
    else
        type(end) = 'P';
    end
else
    index(2 : pP : end) = index(2 : pP : end) + pP;    
end