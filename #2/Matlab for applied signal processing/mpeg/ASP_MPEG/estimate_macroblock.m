function [bcomp, mvv, mvh, errorf] = estimate_macroblock(b, imr, pos, bs, sa, algorithm)



if(nargin == 5)
    algorithm = 'fullsearch';
end


fil = (pos(1)-1)*bs(1)+1;
col = (pos(2)-1)*bs(2)+1;

errorf = 15000*ones(sa(1)*2+1,sa(2)*2+1);

switch(algorithm)
case 'fullsearch'

    b = double(b);
    bb = double(imr(fil:fil+bs(1)-1,col:col+bs(2)-1));
    dmin = sum(sum((b-bb).^2));
    dmin = sum(sum(abs(b-bb)));
    mvv = 0; mvh = 0;
    bcomp = bb;
    
    for r = fil-sa(1):fil+sa(1),
        for c = col-sa(2):col+sa(2),
        
            % block in that position
            bb = double(imr(r:r+bs(1)-1,c:c+bs(2)-1));
        
            % error
            %d = sum(sum((b-bb).^2));
            d = sum(sum(abs(b-bb)));

            errorf(r-fil+sa(1)+1,c-col+sa(2)+1) = d;

            if (d <= dmin),
                dmin = d;
                mvv = r-fil;
                mvh = c-col;
                bcomp = bb;
            end;    
        
        end;
    end;

case 'nstep'
    if(sa ~= [15 15])
        error('Sorry: for now only maximum displacement of 15 allowed...');
    end
	
    is = size(imr);
    
    b = double(b);
    imr = double(imr);
    
    bli = (pos(1) - 1) * bs(1) + 1 : min(is(1), pos(1) * bs(1));
    bco = (pos(2) - 1) * bs(2) + 1 : min(is(2), pos(2) * bs(2));

    bcomp = imr(bli, bco);
    emin = sum(sum(abs(b - imr(bli, bco))));
    mx = 0;mvv=0;
    my = 0;mvh=0;
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
                if(any(brli < 1) | any(brli > is(1)) | any(brco < 1) | any(brco > is(2)))
                    continue;
                end
                e = sum(sum(abs(b - imr(brli, brco))));
                
                % keep error
                errorf(x+sa(1)+1,y+sa(2)+1) = e;
                
                if(emin > e)
                    emin = e;
                    mx = x;
                    my = y;                    
                
                    mvv = mx;
                    mvh = my;
                    bcomp = imr(brli, brco);

                end
            end
        end
        step = floor(step / 2);
    end

otherwise
    error('unknown algorithm');
end

