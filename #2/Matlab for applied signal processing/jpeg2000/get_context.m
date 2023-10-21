function c = get_context(i,j,sigstate,refstate,sbtype)
% INPUTS : 
% - SIGSTATE : matrix with the same size as the processed code-block and
% containing the significant status of each coefficient.
% - REFSTATE : matrix with the same size as the processed code-block and
% containing the "first refinement" status of each coefficient.
% - I,J : row and column indices of the processed coefficient
% - sbtype : kind of subband to which the code-block belongs (LL, HL, LH or
% HH).
%
% OUTPUTS : the context number C according to the following context list : 
% 
% Context list for not yet significant coefficients
% LL and LH subbands     HL subband          HH subband           context
% (vertical high-pass)   (hor. high-pass)    (diag. high-pass)    number
% S(Hi) S(Vi)  S(Di)    S(Hi) S(Vi)  S(Di)   S(Hi+Vi)   S(Di) 
%  2     x      x         x     2      x        x        >=3         9 
%  1     >=1    x         >=1   1      x        >=1      2           8 
%  1     0      >=1       0     1      >=1      0        2           7 
%  1     0      0         0     1      0        >=2      1           6 
%  0     2      x         2     0      x        1        1           5 
%  0     1      x         1     0      x        0        1           4 
%  0     0      >=2       0     0      >=2      >=2      0           3 
%  0     0      1         0     0      1        1        0           2 
%  0     0      0         0     0      0        0        0           1 
%
% Context list for already significant coefficients
% S(Hi)+S(Vi)+S(Di)  First refinement for this coefficient ?       context number 
%      x                            false                            12 
%      >=1                          true                             11 
%      0                            true                             10 
%
% Note : 
% - "x" means "don't care".
% - signs are not entropy coded and are included "as is" in the bit-stream,
%   when the corresponding coefficient becomes significant.

sigborder = zeros(size(sigstate,1)+2,size(sigstate,2)+2);
sigborder(2:end-1,2:end-1)=sigstate;

SHi=sigborder(j+1,i)+sigborder(j+1,i+2);
SVi=sigborder(j,i+1)+sigborder(j+2,i+1);
SDi=sigborder(j,i)+sigborder(j,i+2)+sigborder(j+2,i)+sigborder(j+2,i+2);

if sigstate(j,i)
    if refstate(j,i)
        c=12;
    elseif SHi+SVi+SDi
        c=11;
    else
        c=10;
    end
else
    if sbtype<4
        if sbtype==2
            tmp=SVi;
            SVi=SHi;
            SHi=tmp;
        end
        if SHi==2
            c=9;
        elseif SHi==1
            if SVi>=1
                c=8;
            elseif SDi>=1
                c=7;
            else
                c=6;
            end
        else
            if SVi==2
                c=5;
            elseif SVi==1
                c=4;
            elseif SDi>=2
                c=3;
            elseif SDi==1
                c=2;
            else
                c=1;
            end
        end
    else
        if SDi>=3
            c=9;
        elseif SDi==2
            if SHi+SVi>=1
                c=8;
            else
                c=7;
            end
        elseif SDi==1
            if SHi+SVi>=2
                c=6;
            elseif SHi+SVi==1
                c=5;
            else
                c=4;
            end
        else
            if SHi+SVi>=2
                c=3;
            elseif SHi+SVi==1
                c=2;
            else
                c=1;
            end
        end
    end
end
        

return;