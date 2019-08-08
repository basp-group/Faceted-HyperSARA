%% Fourier reduction parameters
FT2 = @(x) fftshift(fft2(ifftshift(x)));
IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

param_fouRed.enable_klargestpercent = enable_klargestpercent; % 1
param_fouRed.klargestpercent = klargestpercent;
param_fouRed.enable_estimatethreshold = ~enable_klargestpercent;
param_fouRed.gamma = 3;             % By using threshold estimation, the optimal theshold reads as gamma * sigma / ||x||_2
param_fouRed.diagthresholdepsilon = 1e-10; 
param_fouRed.covmatfileexists = 0;
param_fouRed.covmatfile = 'covariancemat.mat';
param_fouRed.fastCov = 1;

fprintf('\nDimensionality reduction...\n');

for i = 1:length(ch)
    if usingBlocking
        for j = 1: length(G{i})
            H{i}{j} = G{i}{j}'*G{i}{j};
            if param_fouRed.enable_estimatethreshold
                param_fouRed.x2 = norm(x0(:, :, i));
                % param_fouRed.dirty2 = norm(Phi_t{i}(y{i})) / sqrt(numel(x0(:,:,i)));
                param_fouRed.dirty2 = norm(operatorPhit(y{i}{j}, G{i}{j}', At),'fro'); % / sqrt(numel(x0(:,:,i)))
                if numel(sigma_noise_ch{i}{j}) == 1
                    param_fouRed.sigma_noise{i}{j} = sigma_noise_ch{i}{j};
                else 
                    tmp = sigma_noise_ch{i}{j};
                    param_fouRed.sigma_noise{i}{j} = mean(tmp(:));
                end
            end
            % no longer invoke function fourierReduction to reduce lambda functions
            % fast way of matrix probing (using psf)
            dirac2D = zeros(Ny, Nx);
            dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;
            PSF = operatorIpsf(dirac2D, A, At, H{i}{j}, [oy*Ny, ox*Nx]);
            covariancemat = FT2(PSF);
            d_mat = abs(real(covariancemat(:)));

            if param_fouRed.enable_klargestpercent
                Wm{i}{j} = (d_mat >= prctile(d_mat,100-param_fouRed.klargestpercent));
            elseif param_fouRed.enable_estimatethreshold
                % Embed the noise
                noise1 = param_fouRed.sigma_noise{i}{j} * (randn(size(G{i}{j}, 1),1) + 1j * randn(size(G{i}{j}, 1), 1));
                rn_tmp = FT2(At(G{i}{j}'*noise1));  % Apply F Phi
                th = param_fouRed.gamma * std(rn_tmp(:)) / param_fouRed.x2;
                % th_dirty = param_fouRed.gamma * std(rn(:)) / param_fouRed.dirty2;
                th_dirty = param_fouRed.gamma * std(rn_tmp(:)) / (param_fouRed.dirty2 / max(d_mat(:)));
                fprintf('\nThe estimate threshold using ground truth is %e \n', th);
                Wm{i}{j} = (d_mat >= th_dirty);
            end
            d_mat = d_mat(Wm{i}{j});
            fprintf('\nThe threshold is %e \n', min(d_mat));

            Ti{i}{j} = max(param_fouRed.diagthresholdepsilon, d_mat);  % This ensures that inverting the values will not explode in computation
            Ti{i}{j} = 1./sqrt(Ti{i}{j});

            yT{i}{j} = dataReduce(y{i}{j}, G{i}{j}', At, Ti{i}{j}, Wm{i}{j});
            rn{i}{j} = dataReduce(noise{i}{j}, G{i}{j}', At, Ti{i}{j}, Wm{i}{j});
            
            G{i}{j} = [];
        end
    
    else
    
        H{i} = G{i}'*G{i};
    
        if param_fouRed.enable_estimatethreshold
            param_fouRed.x2 = norm(x0(:, :, i));
            % param_fouRed.dirty2 = norm(Phi_t{i}(y{i})) / sqrt(numel(x0(:,:,i)));
            param_fouRed.dirty2 = norm(operatorPhit(y{i}, G{i}', At), 'fro'); %/ sqrt(numel(x0(:,:,i))));
            if numel(sigma_noise_ch{i}) == 1
                param_fouRed.sigma_noise = sigma_noise_ch{i};
            else 
                tmp = sigma_noise_ch{i};
                param_fouRed.sigma_noise = mean(tmp(:));
            end
        end

        % no longer invoke function fourierReduction to reduce lambda functions
        % fast way of matrix probing (using psf)
        dirac2D = zeros(Ny, Nx);
        dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;
        PSF = operatorIpsf(dirac2D, A, At, H{i}, [oy*Ny, ox*Nx]);
        covariancemat = FT2(PSF);
        d_mat = abs(real(covariancemat(:)));

        if param_fouRed.enable_klargestpercent
            Mask{i} = (d_mat >= prctile(d_mat,100-param_fouRed.klargestpercent));
        elseif param_fouRed.enable_estimatethreshold
            % Embed the noise
            noise1 = param_fouRed.sigma_noise * (randn(size(G{i}, 1),1) + 1j * randn(size(G{i}, 1), 1));
            rn_tmp = FT2(At(G{i}'*noise1));  % Apply F Phi
            th = param_fouRed.gamma * std(rn_tmp(:)) / param_fouRed.x2;
            % th_dirty = param_fouRed.gamma * std(rn(:)) / param_fouRed.dirty2;
            th_dirty = param_fouRed.gamma * std(rn_tmp(:)) / (param_fouRed.dirty2 / max(d_mat(:)));
            fprintf('\nThe estimate threshold using ground truth is %e \n', th);
            Mask{i} = (d_mat >= th_dirty);
        end
        d_mat = d_mat(Mask{i});
        fprintf('\nThe threshold is %e \n', min(d_mat));

        Sigma{i} = max(param_fouRed.diagthresholdepsilon, d_mat);  % This ensures that inverting the values will not explode in computation
        Sigma{i} = 1./sqrt(Sigma{i});
    
    
%     B{i} = @(x) operatorRPhi(x, A, At, H, W{i}, Sigma{i}, Mask{i}, [Ny, Nx]);
%     Bt{i} = @(x) operatorRPhit(x, A, At, H, W{i}, Sigma{i}, Mask{i}, [Ny, Nx]);
%     FIpsf{i} = @(x) serialise(FT2(operatorIpsf(x, A, At, H{i}, [oy*Ny, ox*Nx], W{i})));  % F * Ipsf, image -> vect
%     FIpsf_t{i} = @(x) operatorIpsf(IFT2(reshape(x, [Ny, Nx])), A, At, H{i}, [Ny, Nx], W{i});  % Ipsf * F^T, vect -> image
   
        yMat = dataReduce(y{i}, G{i}', At, Sigma{i}, Mask{i});
        noiseMat = dataReduce(noise{i}, G{i}', At, Sigma{i}, Mask{i});
        G{i} = [];
        
        if usingReductionPar
            [yT{i}, rn{i}, Ti{i}, Wm{i}] = util_gen_sing_block_structure(yMat, noiseMat, Sigma{i}, Mask{i}, param_sing_block_structure);
        else
            Ti{i} = {Sigma{i}};
            Wm{i} = {Mask{i}};
            yT{i} = {yMat};
            rn{i} = {noiseMat};
        end
    end
    
    R = length(Wm{i});
    aW{i} = cell(R, 1);
    for q = 1:R
        aW{i}{q} = 1./Ti{i}{q};
        epsilont{i}{q} = norm(rn{i}{q});
        epsilons_t{i}{q} = 1.001*epsilont{i}{q}; % data fidelity error * 1.001
    end    
end

clear covariancemat d_mat dirac2D PSF rn y yMat noise noiseMat G Sigma Mask;

fprintf('\nDimensionality reduction is finished\n');
