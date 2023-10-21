function waveshow(sig, W, L)
% Display the 1D DWT as produced by wavedec()
    
    N = length(sig);
    J = log2(N);
    
    %% Drawing initial signal
    n_plot = length(L);
    subplot(n_plot,1,1);
    plot(1:N, sig);
    h = title(sprintf('Initial Signal : N=%i samples, log_2(N)=J=%i', N, J));
    set(h, 'fontsize', 14);
    axis tight;
    
    %% Drawing Wavelet details coefficients
    cW = W;
    for n = (n_plot-1):-1:2;
        j = J - (n_plot - n);
        x = 1:2^(n_plot - n):N;
        
        cwav = cW(end-L(n)+1:end);
        cW = cW(1:end-L(n));
               
        subplot(n_plot,1,(n_plot - n + 1));
        bar(x, cwav);
        axis tight;
        set(gca,'visible','off');
        h = text(x(end) + 10, max(cwav)/2, sprintf('d_%i', j));
        set(h, 'fontsize', 14);
    end
    
    %% Drawing approximation coefficients
    capp = cW(end-L(1)+1:end);
    subplot(n_plot,1,n_plot);
    bar(x, capp);
    axis tight;
    set(gca,'visible','off');
    h = text(x(end) + 10, max(cwav)/2, sprintf('a_%i', j));
    set(h, 'fontsize', 14);
    