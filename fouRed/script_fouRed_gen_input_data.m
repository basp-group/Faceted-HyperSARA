%% generate the sampling pattern
    
if gen_data
    sampling_pattern = 'gaussian+large-holes';

    util_gen_sampling_pattern_config; % Set all parameters
    sparam.N = N; % number of pixels in the image
    sparam.Nox = ox*Nx; % number of pixels in the image
    sparam.Noy = oy*Ny; % number of pixels in the image
    sparam.p = ceil(visibSize/N);

    [~, ~, uw, vw, ~] = util_gen_sampling_pattern(sampling_pattern, sparam);


    if use_symmetric_fourier_sampling
        uw = [uw; -uw];
        vw = [vw; -vw];
    end
    save(uvfile,'uw','vw')
else
    load(uvfile)
end

%% Block structure
if ~usingReductionPar
    param_precond.N = N; % number of pixels in the image
    param_precond.Nox = ox*Nx; % number of pixels in the image
    param_precond.Noy = oy*Ny; % number of pixels in the image
    [aWw] = util_gen_preconditioning_matrix(uw, vw, param_precond);

    %% set the blocks structure
    nWw = ones(length(uw), 1);
    [u, v, ~, uvidx, aW, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);

    %% measurement operator initialization 
    fprintf('Initializing the NUFFT operator\n\n');
    tstart = tic;
    [A, At, G, W, Gw, As, Ats, S] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW, param_nufft);
    tend = toc(tstart);
    fprintf('Initialization runtime: %ds\n\n', ceil(tend));
    R = length(v);
else
    %% If not classical block structure, generate non-block structured data
    %% measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    tstart = tic;
    [A, At, Gw, scale] = op_nufft([vw uw], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], 0);
    tend = toc(tstart);
    fprintf('Initialization runtime: %ds\n\n', ceil(tend));
    
    % use the absolute values to speed up the search
    Gw_a = abs(Gw);
    
    b_l = length(uw);
    % check if eack line is entirely zero
    W = Gw_a' * ones(b_l, 1) ~= 0;
    
    % store only what we need from G
    G = Gw(:, W);    
end

%% sparsity operator definition

[Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
[Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);

%% generate noisy input data in terms of block structure
use_different_per_block_input_snr = 0;
per_block_input_snr_delta = 0;

param_image_var.use_image_variation = 0;

for k = 1:num_tests
    if usingReductionPar
        [y0{k}{1}, y{k}{1}, aY{k}{1}, sigma_noise, noise{k}] = util_gen_input_data_noblock(im, G, W, A, input_snr);
    else
        [y0{k}, y{k}, y0f{k}, yf{k}, aY{k}, input_snr_v{k}, sigma_noise, noise{k}] = util_gen_input_data(im, G, W, A, input_snr, ...
            use_different_per_block_input_snr, per_block_input_snr_delta, uvidx);

        if use_symmetric_fourier_sampling
            y0f{k} = [y0f{k}(uvidx{k}(1:end/2)); conj(y0f{k}(uvidx{k}(1:end/2)))];
            yf{k} = [yf{k}(uvidx{k}(1:end/2)); conj(yf{k}(uvidx{k}(1:end/2)))];
            for j = 1:R
                y{k}{j} = [y{k}{j}(uvidx{k}(1:end/2)); conj(y{k}{j}(uvidx{k}(1:end/2)))];
                y0{k}{j} = [y0{k}{j}(uvidx{k}(1:end/2)); conj(y0{k}{j}(uvidx{k}(1:end/2)))];
                aY{k}{j} = [aY{k}{j}(uvidx{k}(1:end/2)); conj(aY{k}{j}(uvidx{k}(1:end/2)))];
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For dimensionality reduction

if usingReduction
    
    % Continuous-to-gridding projection operator FPhi, singular matrix Sigma, mask matrix (to reduce the dimension)
    [FPhi, FPhit, Sigma, Mask] = fourierReduction(Gw, A, At, [Ny, Nx], param_fouRed);
    % New measurement operator C, new reduced measurement operator B
    [C, Ct, B, Bt] = oper_fourierReduction(FPhi, FPhit, Sigma, Mask, Gw, A, At, [Ny, Nx]);

    evl = op_norm(B, Bt, [Ny, Nx], 1e-4, 200, verbosity);
    
    % Embed the y using the same reduction
    for k = 1:num_tests
        ry = FPhi(y{k}{1});
        yTmat = Sigma.*ry(Mask);
    
        if usingReductionPar
            [yT{k}, T, W] = util_gen_sing_block_structure(yTmat, Sigma, Mask, param_sing_block_structure);
        else
            T = mat2cell([1], 1);
            W = mat2cell(true(size(yTmat)), length(yTmat));
            yT{k} = mat2cell(yTmat, length(yTmat));
        end
    end

    %Bound for the L2 norm
    fprintf('Computing epsilon bound... ');
    tstart1=tic;      
   
    for k = 1:num_tests
        rn = FPhi(noise{k});
        if usingReductionPar        
            for i = 1:length(T)
                epsilonT{k}{i} = norm(T{i} .* rn(W{i}));
                epsilonTs{k}{i} = 1.001*epsilonT{1}{i};
            end
            epsilon{k} = norm(cell2mat(epsilonT{k}));
            epsilons{k} = 1.001*epsilon{k};     % data fidelity error * 1.001
        else
            epsilon{k} = norm(Sigma .* rn(Mask));
            % epsilon = step_epsilon; % set epsilon value BEFORE running this script
            epsilons{k} = 1.001*epsilon{k};     % data fidelity error * 1.001
    %         histogrampeakiness = mean(d12rnnorms)/std(d12rnnorms);
            epsilonT{k} = epsilon;
            epsilonTs{k} = epsilons;
        end
    end
        %%%%%%%%%%%%%%%
    fprintf('Done\n');
    tend1=toc(tstart1);
    fprintf('Time: %e\n', tend1);
else
    T = G;
    Tw = Gw;
    yT = y;
    
    evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-6, 200, verbosity);
    B = A;
    Bt = At;
    
    %Bound for the L2 norm
    fprintf('Computing epsilon bound... ');
    tstart1=tic;  
    for k = 1:num_tests
        [epsilonT{k}, epsilonTs{k}, epsilon{k}, epsilons{k}] = util_gen_L2_bounds(y{k}, ...
        input_snr, [], l2_ball_definition, stopping_criterion, use_same_stop_criterion, param_l2_ball);
    end
    fprintf('Done\n');
    tend1=toc(tstart1);
    fprintf('Time: %e\n', tend1);
end
