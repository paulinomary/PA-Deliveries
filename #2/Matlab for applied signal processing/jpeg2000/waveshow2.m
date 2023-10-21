function waveshow2(W, S, wname)
% Display the 2D DWT as produced by wavedec2()
% Each H,V,D are renormalized together for a nice display
%
% INPUTS :
% * W : vector as produced by wavedec2()
% * S : the structure of coeffs sizes as produced by wavedec2()
% * wname : Wavelet filter name
% 
    
    N = size(S,1) - 2;
    W_img = appcoef2(W, S, wname, N); 
    W_img = W_img ./ max(abs(W_img(:)));
    
    for n = N:-1:1,
        [H,V,D] = detcoef2('all', W, S, n);
        HVD = [H(:), V(:), D(:)];
        mHVD = max(abs(HVD(:)));
        H = H/mHVD;
        V = V/mHVD;
        D = D/mHVD;
        W_img = [W_img, H; V, D];
    end
    
    clf;
    imagesc(W_img);
    axis equal;
    axis tight;
    axis off
    colormap('gray');
    
    hold on;
    for n = (N+1):-1:2,
        X1 = [0 S(n,2)*2];
        Y1 = [S(n,1) S(n,1)];
        
        X2 = [S(n,2) S(n,2)];
        Y2 = [0 S(n,1)*2];
        
        h1 = line(X1, Y1);
        h2 = line(X2, Y2);
        
        set(h1, 'color', 'w');
        set(h2, 'color', 'w');
    end
    
    h = title(sprintf('(Renormalized) %s DWT : %i resolutions', wname, N));
    set(h, 'fontsize', 14);
    colorbar;