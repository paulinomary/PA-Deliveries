function plot_dct_reconstructed(d,vals,type)

% plot_dct_reconstructed - plots recontructed dct matrix
%    plot_dct_reconstructed(d,vals,type) reconstructs dct matrix d zeroing coefficients in vals
%    a new figure is created

switch lower(type),
    case 'zigzag'
        v = zigzag8(d);
        n = 1;
        for N=vals,
            ve = [v(1:N),zeros(1,64-N)];
            dr = izigzag8(ve);
            br = uint8(idct2(dr)+128);
            %subplot(2,3,n);
            imshow(br);
            title(['N = ' int2str(N)]);
            n = n + 1;
        end;
    case 'highest'
        d2 = d(:);
        [v,i] = sort(abs(d2),1,'descend');
        v = d2(i);
        ii(i) = 1:64;
        dct_energy_distribution = cumsum(v.*v)/sum(v.*v);
        n = 1;
        figure;
        for N=vals,
            ve = [v(1:N);zeros(64-N,1)];
            veu = ve(ii);
            dr2 = reshape(veu,8,8);
            br = uint8(idct2(dr2)+128);
            subplot(2,3,n);
            imshow(br);
            title(['N = ' int2str(N) ' (~' int2str(round(dct_energy_distribution(N)*100)) '%)']);
            n = n + 1;
        end;
end


