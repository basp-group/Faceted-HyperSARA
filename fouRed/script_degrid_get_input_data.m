%% generate the sampling pattern
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

param_precond.N = N; % number of pixels in the image
param_precond.Nox = ox*Nx; % number of pixels in the image
param_precond.Noy = oy*Ny; % number of pixels in the image
[aWw] = util_gen_preconditioning_matrix(uw, vw, param_precond);

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
    param_fouRed.dirty2 = norm(Phi_t(y{1}{1}));
    if normalize_data
        param_fouRed.sigma_noise = 1;
    else
        param_fouRed.sigma_noise = sigma_noise;
    end
end

fprintf('\nDimensionality reduction...');

% Reduction using gridding kernel
H = Gw' * Gw;
d = full(abs(diag(H)));

if param_fouRed.enable_klargestpercent
    Mask = (d >= prctile(d,100-klargestpercent));
elseif param.enable_estimatethreshold
    % Embed the noise
    noise = sigma_noise/sqrt(2) * (randn(size(Gw, 1),1) + 1j * randn(size(Gw, 1), 1));
    rn = FT2(At(Gw'*noise));  % Apply F Phi^T
    stdn = std(rn(:));
    param.stdn = stdn;
    th = gamma * stdn / x2;
%     th_dirty = gamma * std(rn(:)) / dirty2;
%     param.maxd = max(d);
%     param.th = th;
%     param.th_dirty = th_dirty;
%     param.singInit = d;
    fprintf('\nThe estimate threshold using ground truth is %e \n', th);
%     fprintf('\nThe estimate threshold using dirty image is %e \n', th_dirty);
    Mask = (d >= th);
    th_per = sum(Mask) / numel(Mask) * 100;
    fprintf('%f%% data are kept\n', th_per);
end
d = d(Mask);
fprintf('\nThe threshold is %e \n', min(d));

d = max(1.e-10, d);  % This ensures that inverting the values will not explode in computation
Sigma = 1./sqrt(d);

FIpsf = @(x) H * A(x);
FIpsf_t = @(x) At(H*x);
fprintf('\nDimensionality reduction is finished');

% Embed the y and noise using the same reduction
for k = 1:num_tests
%     norm_y = norm(y{k}{1})/sqrt(numel(y{k}{1}));
    y_tmp = Gw' * y{k}{1};
    yTmat = Sigma .* y_tmp(Mask);
%     norm_yT = norm(yTmat)/sqrt(numel(yTmat));
%     precond_factor = norm_yT / norm_y;
%     fprintf('Precondition factor=%f\n', precond_factor);
    noise_tmp = Gw' * noise{k}{1};
    noiseMat = Sigma .* noise_tmp(Mask);
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
    aWw_tmp = abs(Gw' * aWw);
    aW = {Sigma .* aWw_tmp(Mask)};
%     aW{q} = 1./T{q};
end

if usingPrecondition
    evl = op_norm(@(x) sqrt(cell2mat(aW)) .* operatorGtPhi(x, H, A, Sigma, Mask), ... 
    @(x) operatorGtPhi_t(sqrt(cell2mat(aW)) .* x, H, At, Sigma, Mask, [oy*Ny, ox*Nx]), ...
        [Ny, Nx], 1e-10, 200, verbosity);
else
    evl = op_norm(@(x) operatorGtPhi(x, H, A, Sigma, Mask), ...
        @(x) operatorGtPhi_t(x, H, At, Sigma, Mask, [oy*Ny, ox*Nx]),...
        [Ny, Nx], 1e-10, 200, verbosity); 
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
