%% generate the sampling pattern
FT2 = @(x) fftshift(fft2(ifftshift(x)));
IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

if gen_data == 0
    fprintf('Loading data from disk ... \n\n');
    load(uvfile)
    
elseif gen_data == 1
    [im, N, Ny, Nx] = util_read_image(image_file_name);

    param_sampling.N = N; % number of pixels in the image
    param_sampling.Nox = ox*Nx; % number of pixels in the image
    param_sampling.Noy = oy*Ny; % number of pixels in the image
%     param_sampling.sigma = sigma_gaussian;
    util_gen_sampling_pattern_config; % Set all parameters

    [~, ~, uw, vw, ~] = util_gen_sampling_pattern(sampling_pattern, param_sampling);


    if use_symmetric_fourier_sampling
        uw = [uw; -uw];
        vw = [vw; -vw];
    end
    
    if save_data_on_disk == 1
        fprintf('Saving new data ... \n');
        
        if save_data_on_disk
            save(uvfile,'uw','vw')
        end
    end
end

%% Plots for test
figure()
plot(uw, vw, '.')

%% measurement operator initialization
fprintf('Initializing the NUFFT operator\n\n');
tstart = tic;
[A, At, Gw, scale] = op_nufft([vw uw], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], 0);
tend = toc(tstart);
fprintf('Initialization runtime: %ds\n\n', ceil(tend));

% use the absolute values to speed up the search
% Gw_a = abs(Gw);

% b_l = length(uw);
% check if eack line is entirely zero
% W = Gw_a' * ones(b_l, 1) ~= 0;

% store only what we need from G
% G = Gw(:, W);    
% end

%% sparsity operator definition

[Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
[Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);

%% generate noisy input data
 
for k = 1:num_tests
    if ~exist('W', 'var')
        W = true(size(Gw, 2), 1);
    end
    [y0{k}{1}, y{k}{1}, ~, sigma_noise, noise{k}{1}] = util_gen_input_data_noblock(im, Gw, W, A, input_snr);       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For natural weighting

if normalize_data
    if numel(sigma_noise) == 1
        Gw = 1./sigma_noise * Gw;      % Whitening G matrix (embed natural weighting in the measurement operator). In reality, this should be done by natural weighting!
        for k = 1:num_tests
            y{k}{1} = 1./sigma_noise * y{k}{1};
            noise{k}{1} = 1./sigma_noise * noise{k}{1};
        end
    else
        if numel(sigma_noise) == length(uw)
            Gw = diag(1./sigma_noise) * Gw;
            for k = 1:num_tests
                y{k}{1} = diag(1./sigma_noise) * y{k}{1};
                noise{k}{1} = diag(1./sigma_noise) * noise{k}{1};
            end
        else
            error('Dimension of natural weights does not match')
        end
    end
end

Phi = @(x) operatorPhi(x, Gw, A);
Phi_t = @(x) operatorPhit(x, Gw', At);

%% For dimensionality reduction
if param_fouRed.enable_estimatethreshold
    param_fouRed.x2 = norm(im);
%     dir2=zeros(Ny,Nx);
%     dir2(ceil((Ny+1)/2),ceil((Nx+1)/2))=1;
%     psf = At(Gw' * Gw * A(dir2));
%     param_fouRed.factor1 = max(abs(psf(:)));
%     param_fouRed.factor2 = sum(abs(psf(:)))/sqrt(N);
%     param_fouRed.factor3 = norm((abs(psf)))/sqrt(N);
%     param_fouRed.dirty2 = norm(At(Gw' * y{1}{1}))/sqrt(N);
    param_fouRed.dirty2 = norm(Phi_t(y{1}{1}))/sqrt(N);
    param_fouRed.sigma_noise = sigma_noise;
end

fprintf('\nDimensionality reduction...');
% psf operator Ipsf, singular value matrix Sigma, mask matrix (to reduce the dimension)
[Ipsf, Mask, Sigma, FIpsf, FIpsf_t] = fourierReduction(Gw, A, At, [Ny, Nx], W, param_fouRed);
% New measurement operator C, new reduced measurement operator B
% [C, Ct, B, Bt] = oper_fourierReduction(Ipsf, Sigma, Mask, [Ny, Nx]);
fprintf('\nDimensionality reduction is finished');

% Embed the y and noise using the same reduction
for k = 1:num_tests
    yTmat = operatorR(y{k}{1}, Phi_t, Sigma, Mask);
    noiseMat = operatorR(noise{k}{1}, Phi_t, Sigma, Mask);
    % Data splitting 
    if usingReductionPar
        [yT{k}, rn{k}, T, W] = util_gen_sing_block_structure(yTmat, noiseMat, Sigma, Mask, param_sing_block_structure);
    else
        T = {Sigma};
        W = {Mask};
        yT{k} = {yTmat};
        rn{k} = {noiseMat};
    end
end

clear noise yTmat noiseMat;

R = length(W);
for q = 1:R
    aW{q} = 1./T{q};
end

if usingPrecondition
    evl = op_norm(@(x) sqrt(cell2mat(aW)) .* operatorRPhi(x, Ipsf, Sigma, Mask, [Ny, Nx]), ...
        @(x) operatorRPhit(sqrt(cell2mat(aW)) .* x, Ipsf, Sigma, Mask, [Ny, Nx]), ...
        [Ny, Nx], 1e-6, 200, verbosity);
else
    evl = op_norm(@(x) operatorRPhi(x, Ipsf, Sigma, Mask, [Ny, Nx]), ...
        @(x) operatorRPhit(x, Ipsf, Sigma, Mask, [Ny, Nx]), [Ny, Nx], 1e-4, 200, verbosity); 
end

%Bound for the L2 norm
fprintf('Computing epsilon bound... ');
tstart1=tic;      

% Embed the noise
for k = 1:num_tests
    % factorized by singular values and compute l2 ball       
    for q = 1:length(T)
        epsilonT{k}{q} = norm(rn{k}{q});
        epsilonTs{k}{q} = 1.001*epsilonT{k}{q};
    end
    epsilon{k} = norm(cell2mat(epsilonT{k}));
    epsilons{k} = 1.001*epsilon{k};     % data fidelity error * 1.001
end

    %%%%%%%%%%%%%%%
fprintf('Done\n');
tend1=toc(tstart1);
fprintf('Time: %e\n', tend1);
