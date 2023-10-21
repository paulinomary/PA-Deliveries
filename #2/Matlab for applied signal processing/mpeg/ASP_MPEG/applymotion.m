function im = applymotion(imr, vv, vh, bs)

%applymotion   Application of motion vectors to a reference image.
%
%   im = applymotion(imr, vv, vh, bs) applies motion vectors vv and hh to a 
%   reference image, thus obtaining a prediction of the current image.  Returns
%   the prediction 'im' of the current image obtained by applying to the
%   reference image 'imr' motion vectors whose coordinates are given in 'vv'
%   (vertical component) and 'vh' (horizontal component).  The motion vectors
%   are assumed to apply to each block of the image, where the block size is
%   given by vector 'bs' (typicaly [16 16]).
%
%   The motion vectors must not refer to the outside of the reference image. 
%   They normally have integer coordinates.  If the coordinates are not
%   integer, bilinear interpolation is performed.  'vv' and 'vh' must be
%   matrices of the same size.
%
%   See also estimatemotion, applymotionb.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if nargin ~= 4
    error('wrong number of arguments');
end

is = size(imr);

if ~ismatrix(vv) | ~ismatrix(vh)
    error('vv and vh must be matrices');
elseif ~isrow(bs) | length(bs) ~= 2 | ~myisinteger(bs)
    error('bs must be a length 2 integer row vector');
elseif ~samesize(vv, vh)
    error('vv and vh must have the same size');
end

vs = ceil(is ./ bs(1:2));

if ~equal(vs, size(vv))
    error('bs incompatible with im and vv/vh sizes.');
end
    
im = imr;

if ~isinteger(vv) | ~isinteger(vh)
    imr = double(imr);
    imr(end + 1, end + 1) = 0;
end

for li = 1 : vs(1)
    for co = 1 : vs(2)
        bli = (li - 1) * bs(1) + 1 : min(is(1), li * bs(1));
        bco = (co - 1) * bs(2) + 1 : min(is(2), co * bs(2));
            
        if ~isinteger(vv) | ~isinteger(vh)
            brli = bli + floor(vv(li, co));
            brco = bco + floor(vh(li, co));
	
            if any(brli <= 0) | any(is(1) < brli) | ...
                any(brco <= 0) | any(is(2) < brco)
                error('motion vectors must refer to inside of reference image');
            end

            brlif = vv(li, co) - floor(vv(li, co));
            brcof = vh(li, co) - floor(vh(li, co));
        
            if isa(im, 'uint8')
                im(bli, bco, :) = round((1 - brcof) * ((1 - brlif) * imr(brli, brco) + brlif * imr(brli + 1, brco)) + ...
                    brcof * ((1 - brlif) * imr(brli, brco + 1) + brlif * imr(brli + 1, brco + 1)));
            else
                im(bli, bco, :) = (1 - brcof) * ((1 - brlif) * imr(brli, brco) + brlif * imr(brli + 1, brco)) + ...
                    brcof * ((1 - brlif) * imr(brli, brco + 1) + brlif * imr(brli + 1, brco + 1));
            end
        else
            brli = bli + vv(li, co);
            brco = bco + vh(li, co);
	
            if any(brli <= 0) | any(is(1) < brli) | ...
                any(brco <= 0) | any(is(2) < brco)
                error('motion vectors must refer to inside of reference image');
            end

            im(bli, bco, :) = imr(brli, brco, :);    
        end
    end
end

