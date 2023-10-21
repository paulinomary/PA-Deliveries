function appshow(sig, W, L, wname)
% Display the approximation coefficients of 1D DWT as produced by wavedec()

N = length(sig);
J = log2(N);
J0 = J - length(L) + 2;

% Drawing initial signal
n_plot = length(L) - 1;
subplot(n_plot,1,1);
plot(1:N, sig);
h = title(sprintf('Initial Signal : N=%i samples, log_2(N)=J=%i', N, J));
set(h, 'fontsize', 14);
axis tight;

% Initialization
cW = W;

% Drawing the coarsest approximation
capp = cW(1:L(1));
cW = cW((L(1)+1):end);
subplot(n_plot,1,n_plot);
x = 1:(N/L(1)):N;
bar(x, capp);
axis tight;
set(gca,'visible','off');
h = text(x(end) + 10, max(capp)/2, sprintf('a_%i', J0));
set(h, 'fontsize', 14);


% Drawing the other approximation coefficients
for n = 2:(n_plot-1);

    cdet = cW(1:L(n));
    cW = cW((L(n)+1):end);
    capp = idwt(capp, cdet, wname);

    j = J0 + n - 1;
    x = 1:(N/(2*L(n))):N;

    subplot(n_plot, 1, n_plot - n + 1);
    bar(x, capp);
    axis tight;
    set(gca,'visible','off');
    h = text(x(end) + 10, max(capp)/2, sprintf('a_{%i}', j));
    set(h, 'fontsize', 14);
end