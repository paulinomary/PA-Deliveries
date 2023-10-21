function [vv, vh, errors] = estimatemotion(im, imr, bs, dmax, algorithm)

%estimatemotion   Block matching motion estimation.
%
%   [vv, vh, errors] = estimatemotion(im, imr, bs, dmax, algorithm])
%   estimates of the motion vectors of blocks obtained using the block
%   matching technique.  Estimates motion vectors with components 'vv'
%   (vertical) and 'vh' (horizontal) for each block with size 'bs' of image
%   'im' with reference to image 'imr'.  The motion vectors will have a
%   maximum displacement of 'dmax' (a vector, since the maximum
%   displacements may be diferent in each  direction).  The algorithm used
%   by default is 'fastfullsearch'.  It may be  changed to 'nstep', where a
%   n-step refination search is performed, or to  'fullsearch'.  Both
%   'nstep' and 'fullsearch' minimize the sum of the absolute differences.
%   'fastfullsearch' minimizes the sum of squares. 
%
%   'errors' is a matrix containing the error for each macroblock.  In the
%   case of 'nstep' and 'fullsearch', the error is the sum of the absolute
%   differences.  Int the case of 'fastfullsearch', it is the sum of
%   squares (notice: not the square root of the sum of squares!).
%
%   'fastfullsearch' is about 13 times faster than 'fullsearch'.  'nstep' is
%   about 2 times faster than 'fastfullsearch'.  (Results otained for CIF
%   images).
%
%   The estimated motion vectors never point outside of the reference image.
%
%   The PSNR between the original and the motion compensated image is usually
%   considerable worse with 'nstep' than with 'fullsearch', which in turn is 
%   usually usually slightly worse than with 'fastfullsearch'.
%
%   Whenever two vectors produce the same smallest error, the vector with
%   the smaller norm is chosen.
%
%   Bugs: for now 'nstep' only works with 'dmax = [15 15]'.
%
%   See also estimatemotion, applymotionb.
%
%Copyright (C) 2001-2002 - Manuel Menezes de Sequeira (ISCTE).

if nargin ~= 4 & nargin ~= 5
    error('wrong number of arguments');
end

is = size(im);

if ~isrow(bs) | length(bs) ~= 2 | ~myisinteger(bs)
    error('bs must be a length 2 integer row vector');
elseif ~isrow(dmax) | length(dmax) ~= 2 | ~myisinteger(dmax)
    error('dmax must be a length 2 integer row vector');
elseif ~samesize(im, imr)
    error('im and imr must have the same size');
end

vs = ceil(is ./ bs(1:2));

if(nargin == 4)
    algorithm = 'fullsearch';
end

is = size(im);

im = double(im);
imr = double(imr);

switch(algorithm)
    
    case 'fullsearch'
	% Build scan order:
	i = -dmax(1) : dmax(1); 
	j = -dmax(2) : dmax(2); 
	[ii, jj] = ndgrid(i, j);
	distances = ii(:) .^ 2 + jj(:) .^ 2; 
	[sd, indices] = sort(distances);
	
	for li = 1 : vs(1)
        for co = 1 : vs(2)
            bli = (li - 1) * bs(1) + 1 : min(is(1), li * bs(1));
            bco = (co - 1) * bs(2) + 1 : min(is(2), co * bs(2));
            emin = sum(sum(abs(im(bli, bco) - imr(bli, bco))));
            mx = 0;
            my = 0;
            for n = 1 : length(indices),
                x = ii(indices(n));
                y = jj(indices(n));
                brli = bli + x;
                brco = bco + y;
                if(any(brli < 1) | any(brli > is(1)) | ...
                   any(brco < 1) | any(brco > is(2)))
                    continue;
                end
                e = sum(sum(abs(im(bli, bco) - imr(brli, brco))));
	
                if(emin > e)
                    emin = e;
                    mx = x;
                    my = y;
                end
            end
            vv(li, co) = mx;
            vh(li, co) = my;
            errors(li, co) = emin;
        end
	end
case 'nstep'
    if(dmax ~= [15 15])
        error('Sorry: for now only maximum displacement of 15 allowed...');
    end
	
	for li = 1 : vs(1)
        for co = 1 : vs(2)
            bli = (li - 1) * bs(1) + 1 : min(is(1), li * bs(1));
            bco = (co - 1) * bs(2) + 1 : min(is(2), co * bs(2));
            emin = sum(sum(abs(im(bli, bco) - imr(bli, bco))));
            mx = 0;
            my = 0;
            step = 8;
            while step ~= 0,
                rx = mx;
                ry = my;
                for i = -step:step:step,
                    for j = -step:step:step,
                        x = rx + i;
                        y = ry + j;
                        brli = bli + x;
                        brco = bco + y;
                        if(any(brli < 1) | any(brli > is(1)) | ...
                           any(brco < 1) | any(brco > is(2)))
                            continue;
                        end
                        e = sum(sum(abs(im(bli, bco) - imr(brli, brco))));
	
                        if(emin > e)
                            emin = e;
                            mx = x;
                            my = y;
                        end
                    end
                end
                step = floor(step / 2);
            end
            vv(li, co) = mx;
            vh(li, co) = my;
            errors(li, co) = emin;
        end
	end
case 'fastfullsearch'

    imr2 = filter2(ones(bs), imr .^ 2);

    for li = 1 : vs(1)
        for co = 1 : vs(2)
            % Coordinates of pixels in current block and current block size (tests 
            % necessary since bottom and right blocks may have to be smaller):
            if li * bs(1) > is(1)
                bli = (li - 1) * bs(1) + 1 : is(1);
                bsli = is(1) - bli(1) + 1;
                special_block = true;
            else
                bli = (li - 1) * bs(1) + 1 : li * bs(1);
                bsli = bs(1);
                special_block = false;
            end
            
            if co * bs(2) > is(2)
                bco = (co - 1) * bs(2) + 1 : is(2);
                bsco = is(2) - bco(1) + 1;
                special_block = true;
            else
                bco = (co - 1) * bs(2) + 1 : co * bs(2);
                bsco = bs(2);
                special_block = false;
            end
            
            % Extract current block:
            block = im(bli, bco);
            
            % Get corresponding superblock on reference image.  This superblock
            % includes all pixels which will be used for comparisons.  Calculate 
            % offsets so that motion vectors can be easily calculated:
            if bli(1) <= dmax(1)
                bliroffset = dmax(1) - bli(1) + 1;
            else
                bliroffset = 0;
            end
            blirfirst = bli(1) - dmax(1) + bliroffset;
            
            if bco(1) <= dmax(2)
                bcoroffset = dmax(2) - bco(1) + 1;
            else
                bcoroffset = 0;
            end
            bcorfirst = bco(1) - dmax(2) + bcoroffset;
            
            blirlast = min(bli(end) + dmax(1), is(1));
            bcorlast = min(bco(end) + dmax(2), is(2));
            
            blir = blirfirst : blirlast;
            bcor = bcorfirst : bcorlast;
            
            % Extract reference block:
            blockr = imr(blir, bcor);
            
            if special_block
                block2r = filter2(ones(size(block)), blockr .^ 2, 'valid');
            else
                % Extract the relevant part of the squared and filtered reference image:
                
                blirvalid = blirfirst + floor((bsli - 1) / 2) : blirlast - ceil((bsli - 1) / 2);
                bcorvalid = bcorfirst + floor((bsco - 1) / 2) : bcorlast - ceil((bsco - 1) / 2);
                
                block2r = imr2(blirvalid, bcorvalid);
            end
            
            % Filter returning only valid values (no zero-padding):
            x = block2r - 2 * filter2(block, blockr, 'valid');
            
            lineorigin = dmax(1) + 1 - bliroffset;
            columnorigin = dmax(2) + 1 - bcoroffset;
            
            i = (1 : size(x, 1)) - lineorigin; 
            j = (1 : size(x, 2)) - columnorigin; 
            [ii, jj] = ndgrid(i, j);
            distances = ii(:) .^ 2 + jj(:) .^ 2; 
            [sd, indices] = sort(distances);
            
            [minimum, n] = min(x(indices));
            
            vv(li, co) = ii(indices(n));
            vh(li, co) = jj(indices(n));
        end
    end
    if nargout > 2
        errors = blkproc((applymotion(imr, vv, vh, bs) - im) .^ 2, bs, 'sum(x(:))');
    end
    
otherwise
    error('unknown algorithm');
end


