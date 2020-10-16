function [yT, H, W, T, aW, Wm, epsilon, epsilons] = Fourier_reduction_ch(y, G, reduction_version, fouRed_gamma, param_fouRed, save_fouRed, load_fouRed)

FT2 = @(x) fftshift(fft2(ifftshift(x))) / sqrt(numel(x));
precision = 1e-16;

if load_fouRed
    DRfilename = ['./data/DR_fouRed', num2str(reduction_version), '_percent', num2str(fouRed_gamma), '.mat'];
    load(DRfilename, 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'epsilon', 'epsilons');
else

    for i = 1:length(ch)

        for j = 1:length(G{i})
            Hl = G{i}{j}' * G{i}{j};
            yw = G{i}{j}' * y{i}{j};   % Grided vis

            G{i}{j} = [];   % save RAM

            % remove small values according to numeric precision
            peak = max(max(abs(Hl)));    
            Hl = Hl .* (abs(Hl) > peak * precision);

            fprintf('Channel: %d, block: %d, initial H matrix memory: \n', i, j)
            whos Hl

            if reduction_version == 1    
            %     fast matrix probing (using psf)
                dirac2D = zeros(Ny, Nx);
                dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;
                PSF = operatorIpsf(dirac2D, A, At, Hl, [oy*Ny, ox*Nx]);
                covariancemat = FT2(PSF);
                d_mat = abs((covariancemat(:)));
            elseif reduction_version == 2
                d_mat = full(abs(diag(Hl)));
            end

            if param_fouRed.enable_klargestpercent
    %             th = 0;
    %             Mask = (d_mat > th);
                Mask = (d_mat >= prctile(d_mat,100-param_fouRed.klargestpercent));
            elseif param_fouRed.enable_estimatethreshold
            %     estimate threshold
                d_mat_sort = sort(d_mat);
                d_mat_sort_cumsum = cumsum(d_mat_sort);
                d_mat_sum = d_mat_sort_cumsum(end); % whole energy of d_mat
                th_energy = d_mat_sum * (1 - prob); % threshold according to the k-sigma rule
                th = d_mat_sort(find(d_mat_sort_cumsum >= th_energy, 1, 'first'));
                Mask = (d_mat >= th);
            end

            ind_nz = d_mat > 0;       % non-zero singular values
            kept = sum(Mask(:));
            total = numel(Mask);
            percentage = kept / total * 100;
            fprintf('\nChannel %d, block %d: %d non-zero singular values, %d over %d, or %f%% of data are kept, threshold=%e\n', i, j, sum(ind_nz), kept, total, percentage, th);

            d_mat = d_mat(Mask);
            Hl = Hl(Mask,:);

            fprintf('Reduced H matrix memory: \n')
            whos Hl

            Tl = d_mat;
            Tl = 1./sqrt(Tl);
            Wml = Mask;

            clear d_mat Mask

            if reduction_version == 1
                if usingPrecondition
                    aWl = Tl;
                else
                    aWl = 1;
                end
                temp = zeros(size(W{i}{j}, 1), 1);
                temp(W{i}{j}) = yw;
                im = FT2(real(At(yw)));
                im = im(:);
                yTl = Tl .* im(Wml);
            elseif reduction_version == 2
                if usingPrecondition
                    aWl = Tl;
                else
                    aWl = 1;
                end
                yTl = Tl.*yw(Wml,:);
            end
        end

        H{i}{j} = Hl;
        yT{i}{j} = yTl;
        T{i}{j} = Tl;
        aW{i}{j} = aWl;
        Wm{i}{j} = Wml;

    end
    
    if save_fouRed
        DRfilename = ['./data/DR_fouRed', num2str(reduction_version), '_percent', num2str(fouRed_gamma), '.mat'];
        if ~isfile(DRfilename)
            save(DRfilename, '-v7.3', 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'epsilon', 'epsilons');
        end
    end
    
end