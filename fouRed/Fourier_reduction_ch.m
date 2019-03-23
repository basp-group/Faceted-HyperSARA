%% Fourier reduction parameters
FT2 = @(x) fftshift(fft2(ifftshift(x)));
IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

param_fouRed.enable_klargestpercent = 1;
param_fouRed.klargestpercent = klargestpercent;
param_fouRed.enable_estimatethreshold = 0;
param_fouRed.gamma = 3;             % By using threshold estimation, the optimal theshold reads as gamma * sigma / ||x||_2
param_fouRed.diagthresholdepsilon = 1e-10; 
param_fouRed.covmatfileexists = 0;
param_fouRed.covmatfile = 'covariancemat.mat';
param_fouRed.fastCov = 1;

fprintf('\nDimensionality reduction...');
for i = 1:length(ch)
    
    Phi{i} = @(x) operatorPhi(x, G{i}, A, [oy*Ny ox*Nx], W{i});
    Phi_t{i} = @(x) operatorPhit(x, G{i}', At, [oy*Ny ox*Nx], W{i});
    
    if param_fouRed.enable_estimatethreshold
        param_fouRed.x2 = norm(x0(:, :, i));
        param_fouRed.dirty2 = norm(Phi_t{i}(y{i})) / sqrt(numel(x0(:,:,i)));
        if numel(sigma_noise{i}) == 1
            param_fouRed.sigma_noise = sigma_noise{i};
        else 
            tmp = sigma_noise{i};
            param_fouRed.sigma_noise = mean(tmp(:));
        end
    end
    
    % psf operator Ipsf, singular value matrix Sigma, mask matrix (to reduce the dimension)
    [Ipsf{i}, Mask{i}, Sigma{i}, FIpsf{i}, FIpsf_t{i}] = fourierReduction(G{i}, A, At, [Ny, Nx], W{i}, param_fouRed);
    
    B{i} = @(x) operatorRPhi(x, Ipsf{i}, Sigma{i}, Mask{i}, [Ny, Nx]);
    Bt{i} = @(x) operatorRPhit(x, Ipsf{i}, Sigma{i}, Mask{i}, [Ny, Nx]);
    
    yMat = operatorR(y{i}, Phi_t{i}, Sigma{i}, Mask{i});
    noiseMat = operatorR(noise{i}, Phi_t{i}, Sigma{i}, Mask{i});
    
    if usingReductionPar
        [yT{i}, rn{i}, T{i}, W{i}] = util_gen_sing_block_structure(yMat, noiseMat, Sigma{i}, Mask{i}, param_sing_block_structure);
    else
        T{i} = {Sigma{i}};
        W{i} = {Mask{i}};
        yT{i} = {yMat};
        rn{i} = {noiseMat};
    end
    
    R = length(W{i});
    aW{i} = cell(R, 1);
    for q = 1:R
        aW{i}{q} = 1./T{i}{q};
        epsilont{i}{q} = norm(rn{i}{q});
        epsilons_t{i}{q} = 1.001*epsilont{i}{q};     % data fidelity error * 1.001
    end
    
end

clear y noise yTmat noiseMat;

fprintf('\nDimensionality reduction is finished');