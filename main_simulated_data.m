function main_simulated_data(image_name, nChannels, Qx, Qy, Qc, p, input_snr, ...
    algo_version, window_type, ncores_data, ind, overlap_size, ...
    flag_generateCube, flag_generateCoverage, flag_generateVisibilities, flag_generateUndersampledCube, ...
    flag_computeOperatorNorm, flag_solveMinimization, ...
    cube_path, coverage_path)
% Main script to run the faceted HyperSARA approach on synthetic data.
% 
% This script generates synthetic data and runs the faceted HyperSARA 
% approach to reconstruct an :math:`N \times L` wideband image cube.
% 
% Args:
%     image_name (string): name of the reference synthetic image (from the 
%     data/ folder)
%     nChannels (int): number of channels
%     Qx (int): number of spatial facets along axis 2 (x)
%     Qy (int): number of spatial facets along axis 1 (y)
%     Qc (int): number of spectral facets
%     p (double): [description]
%     input_snr (double): input SNR value (in dB)
%     algo_version (string): selected version of the solver:
%        - 'standard'      overlap of the faceted nuclear norm equal to the
%                          one needed for the faceted SARA dictionary, w/o 
%                          spatial weights (apodization window)
%        - 'standard2'     alternative parallellization scheme compared to
%                          'standard' (gather image and data on the same 
%                          nodes)
%        - 'cst'           constant overlap for the faceted nuclear norm,  
%                          w/o spatial weights (apodization window)
%        - 'weighted'      overlap of the faceted nuclear norm equal to the 
%                          one needed for the faceted SARA dictionary,  
%                          using spatial weights (apodization window)
%        - 'cst_weighted'  constant overlap taken for the faceted nuclear
%                          norm, using spatial weights (apodization window)
%        - 'no_overlap'    no overlap for the faceted nuclear norm prior
%     window_type (string): type of apodization window considered for the 
%                           faceted nuclear norm prior. Only active with 
%                           the following versions of the algorithm:
%        - 'cst_weighted'
%        - 'weighted'
%     ncores_data (int): number of cores handlig the data fidelity terms 
%                        ("data cores"). The total number of cores is 
%                        Qx*Qy + ncores_data + 1
%     ind (int): index of the spectral facet to be reconstructed (set to -1
%                to deactivate spectral faceting)
%     overlap_size (int): number of overlapping pixels between contiguous 
%                         facets (only active for  the 'cst' and 
%                         'cst_overlap' versions of the solver)
%     flag_generateCube (bool): flag specifying whether the ground truth image 
%                           cube needs to be generated or loaded from an 
%                           existing .fits file
%     flag_generateCoverage (bool): flag specifying whether the uv-coverage 
%     needs to be generated or loaded from an existing .fits file
%     flag_generateVisibilities (bool): flag specifying whether the 
%     visibilities need to be generated or loaded from an existing .mat 
%     file
%     flag_generateUndersampledCube (bool): flag to generate an undersampled 
%     version (by a factor 4) of the ground-truth wideband image cube
%     flag_computeOperatorNorm (bool): compute norm of the measurement 
%     operator, or load it from an existing .mat file (e.g., computed from 
%     a previous run)
%     flag_solveMinimization (bool): flag triggering the faceted HyperSARA 
%     solver
%     cube_path (string): path and name of the wideband cube .fits file 
%     (w/o file extension)
%     coverage_path (string): path and name of the uv-coverage .fits file 
%     (w/o file extension) 

% max_iter
% max_reweight

%% PARAMETERS FOR DEBUGGING

% image_name = 'W28_512_m39';
% nchannels = 60; % (needs to be > 60 to avoid bugs with the version implemented by Abdullah)
% Qx = 2;
% Qy = 1;
% Qc = 1;
% p = 0.3; % percentage
% input_snr = 40; % input SNR (in dB)
% algo_version = 'cst_weighted';
% window_type = 'triangular';
% ncores_data = 1; %number of cores assigned to the data fidelity terms (groups of channels)
% ind = 1;  % index from "spectral" facet
% generate_cube = true;
% generate_coverage = true;
% generate_visibilities = true;
% generate_undersampled_cube = true; % Default 15 channels cube with line emissions
% compute_operator_norm = false;
% solve_minimization = true;
% cube_path = strcat('data/', image_name, '_L', num2str(nchannels));
% coverage_path = 'data/uv_coverage_p=1';
% 
% % if strcmp(algo_version, 'cst_weighted') || strcmp(algo_version, 'cst')
% %     nlevel = 4;
% %     overlap_size = (power(2, nlevel)-1)*(2*8 - 2); % assuming db8 largest wavelet filter
% % end
% 
% overlap_size = 256;
% 
% % % TODO: add warm-restart for this version of the main script

%%
format compact;

disp(['Reference image: ', image_name]);
disp(['nchannels: ', num2str(nChannels)]);
disp(['Number of facets Qy x Qx : ', num2str(Qy), ' x ', num2str(Qx)]);
disp(['Number of spectral facets Qc : ', num2str(Qc)]);
disp(['Overlap size: ', num2str(overlap_size)]);
disp(['Number of data points p per frequency (as a fraction of image size): ', num2str(p)]);
disp(['Input SNR: ', num2str(input_snr)]);
disp(['Generating image cube: ', num2str(flag_generateCube)]);
disp(['Generating coverage: ', num2str(flag_generateCoverage)]);
disp(['Generating visibilities: ', num2str(flag_generateVisibilities)]);

addpath lib/generate_data/
addpath lib/operators/
addpath lib/nufft/
addpath lib/utils/
addpath lib/sdwt2/
addpath data/
addpath src/
addpath src/spmd
addpath src/spmd/weighted
addpath src/spmd/standard
addpath src/spmd/no_overlap

% setting paths to results and reference image cube
coverage_path = strcat(coverage_path, '.fits');
cube_path = strcat(cube_path, '.fits');
data_path = 'data/';
results_path = fullfile('results/', image_name);
reference_cube_path = fullfile(data_path, strcat(image_name, '.fits'));
mkdir(data_path)
mkdir(results_path)

%% 
rng(1234)
% T = 1500; % to be set depending on the value of p
% hrs = 6;
% kernel = 'minmax:tuned'; % 'kaiser', 'minmax:tuned'
generate_eps_nnls = false;
save_data = true; 

%% Generate/load ground-truth image cube
f = linspace(1,2,nChannels);
if flag_generateCube
    % frequency bandwidth from 1 to 2 GHz
    emission_lines = 0; % insert emission lines on top of the continuous spectra
    [x0,X0] = Generate_cube(reference_cube_path,f,emission_lines);
    [Ny, Nx, nChannels] = size(x0);
    if flag_generateUndersampledCube
        % undersampling factor for the channels
        unds = 4; % take 1/unds images
        [x0,X0,f,nChannels] = Generate_undersampled_cube(x0,f,Ny,Nx,nChannels,unds);
    end
    fitswrite(X0.', cube_path);
    fitsdisp(cube_path)
else
    X0 = fitsread(cube_path).';
    Nx = sqrt(size(X0, 1));
    Ny = Nx;
    nChannels = size(X0, 2);
    x0 = reshape(X0, [Ny, Nx, nChannels]);
end

%% Generate spectral facets (interleaved sampling)
id = interleaved_facets(Qc, nChannels);
if ind > 0
    x0 = x0(:,:,id{ind});
    nChannels = size(x0,3);
    f = f(id{ind});
    X0 = reshape(x0,Nx*Ny,nChannels);
end
channels = 1:nChannels;

%% Setup name of results file
data_name_function = @(nchannels) strcat('y_N=',num2str(Nx),'_L=', ...
    num2str(nchannels),'_p=',num2str(p),'_snr=', num2str(input_snr),'.mat');
results_name_function = @(nchannels) strcat('fhs_', algo_version,'_',window_type,'_N=',num2str(Nx), ...
    '_L=',num2str(nchannels),'_p=',num2str(p), ...
    '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc), ...
    '_overlap=', num2str(overlap_size), ...
    '_snr=', num2str(input_snr),'.mat');

data_name = data_name_function(nChannels);
results_name = results_name_function(nChannels);

%% Default config parameters (can be kept as is)
% gridding parameters
N = Nx * Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

% preconditioning parameters
param_precond.N = N;       % number of pixels in the image
param_precond.Nox = ox*Nx; % number of Fourier points (oversampled plane)
param_precond.Noy = oy*Ny;
param_precond.gen_uniform_weight_matrix = 1; % set weighting type
param_precond.uniform_weight_sub_pixels = 1;

% block structure
param_block_structure.use_density_partitioning = 0;
param_block_structure.density_partitioning_no = 1;
param_block_structure.use_uniform_partitioning = 0;
param_block_structure.uniform_partitioning_no = 4;
param_block_structure.use_equal_partitioning = 1;
param_block_structure.equal_partitioning_no = 1;
param_block_structure.use_manual_frequency_partitioning = 0;
% sparam.fpartition = [pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [0, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.25*pi, 0, 0.25*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-64/256*pi, 0, 64/256*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi/4), 0, icdf('norm', 0.75, 0, pi/4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_block_structure.use_manual_partitioning = 0;

%% Generate/load uv-coverage, setup measurement operator
% generating u-v coverage
if flag_generateCoverage
    cov_type = 'vlaa';
    dl = 1.1;
    hrs = 5;
    na = 27; % for vlaa
    M = na*(na-1)/2;
    % Fixing Mt = 0.5 N, take T = 0.5 N / M : na = 27 for vla
    T = floor(p*(Nx*Ny)/M); % should be > 1
    [u, v, ~] = generate_uv_coverage(T, hrs, dl, cov_type);
    % -> in the c++ code, do u.col(f) = u*freq[f]/freq[0]
    % with freq = linspace(1, 2, L)
    % om = [v(:), u(:)];
    u = u(:)*f(1)/f(end);
    v = v(:)*f(1)/f(end);
    fitswrite([u, v, ones(numel(u), 1)], coverage_path)
    fitsdisp(coverage_path);
else
    uvw = fitsread(coverage_path);
    u = uvw(:, 1);
    v = uvw(:, 2);
    % u = u*f(1)/f(end);
    % v = v*f(1)/f(end);
    disp('Coverage loaded successfully')
    clear uvw
end

% setup measurement operator
for i = 1:nChannels
    uw = (f(i)/f(1)) * u;
    vw = (f(i)/f(1)) * v;
    
    % compute uniform weights (sampling density) for the preconditioning
    aWw = util_gen_preconditioning_matrix(uw, vw, param_precond);
    
    % set the weighting matrix, these are the natural weights for real data
    nWw = ones(length(uw), 1);
    
    % set the blocks structure
    [u1, v1, ~, ~, aW{i}, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    [A, At, G{i}, W{i}] = op_p_nufft([v1 u1], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
end

%% Free memory
clear u v u1 v1 uw vw aWw nW nWw param_block_structure param_precond;

%% Generate/load visibilities
if flag_generateVisibilities
    param_l2_ball.type = 'sigma';
    param_l2_ball.sigma_ball = 2;
    [y0, y, Nm, sigma_noise] = util_gen_measurements(x0, G, W, A, input_snr);
    [epsilon,epsilons] = util_gen_data_fidelity_bounds(y, Nm, param_l2_ball, sigma_noise); 

    if save_data
        save(fullfile(data_path,data_name), '-v7.3', 'y0', 'y', 'epsilon', 'epsilons', 'sigma_noise');
    end
else
    load(fullfile(data_path,data_name), 'y', 'epsilons');
    disp('Data loaded successfully')
end

clear y0 Nm;

%% Compute operator norm
if flag_computeOperatorNorm
    % Compute full measurement operator spectral norm
    F = afclean( @(x) HS_forward_operator_precond_G(x, G, W, A, aW));
    Ft = afclean( @(y) HS_adjoint_operator_precond_G(y, G, W, At, aW, Ny, Nx));
    Anorm = pow_method_op(F, Ft, [Ny Nx nChannels]);
    save(fullfile(results_path,strcat('Anorm_N=',num2str(Nx), ...
        '_L=',num2str(nChannels),'_p=',num2str(p),'_snr=', num2str(input_snr), '.mat')),'-v7.3', 'Anorm');
else
    load(fullfile(results_path,strcat('Anorm_N=',num2str(Nx),'_L=', ...
        num2str(nChannels),'_p=',num2str(p),'_snr=', num2str(input_snr), '.mat')));
end
    
%% Generate initial epsilons by performing imaging with NNLS on each data block separately
if generate_eps_nnls
    % param_nnls.im = im; % original image, used to compute the SNR
    param_nnls.verbose = 2; % print log or not
    param_nnls.rel_obj = 1e-3; % stopping criterion
    param_nnls.max_iter = 1000; % max number of iterations
    param_nnls.sol_steps = [inf]; % saves images at the given iterations
    param_nnls.beta = 1;
    % solve nnls per block
    for i = 1:nChannels
        eps_b{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            % printf('solving for band %i\n\n',i)
            [~,eps_b{i}{j}] = fb_nnls_blocks(yb{i}{j}, A, At, G{i}{j}, W{i}{j}, param_nnls);
        end
    end
    mkdir('data/')
    save('data/eps.mat','-v7.3', 'eps_b');
    % load('data/eps.mat');
end

%% Solver
if flag_solveMinimization
    % Definition of the SARA dictionary
    nlevel = 4; % depth of the wavelet decompositions
    %! always specify Dirac basis ('self') in last position if used in the
    %! SARA dictionary 
    wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; 
    L = [2*(1:8)'; 0]; % length of the filters (0 corresponding to the 'self' basis)
    
    % estimate mu: ratio between nuclear and l21 norm priors applied to
    % the dirty image
    dirty_image = zeros([Ny, Nx, nChannels]);
    for l = 1:nChannels
        temp = zeros(oy*Ny*ox*Nx, 1);
        for b = 1:numel(G{l})
            temp(W{l}{b}) = temp(W{l}{b}) + G{l}{b}' * (sqrt(aW{l}{b}) .* y{l}{b});
        end
        dirty_image(:,:,l) = At(temp);
    end
    [~,S0,~] = svd(reshape(dirty_image, [Nx*Ny, nChannels]),'econ');
    nuclear_norm = sum(abs(diag(S0)));
    
    %TODO: do it directly in parallel with faceted SARA
    dwtmode('zpd')
    [~, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    [~, s] = n_wavelet_coefficients(L(1:end-1), [Ny, Nx], 'zpd', nlevel);
    s = s+N; % total number of SARA coefficients
    Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,s);
    %! the number of elements reported below is only valid for the periodic
    %('per') boundary condition (implicitly assumed to be used in HS_forward_sparsity)
    % TODO: enable other different boundary conditions
    l21_norm = sum(sqrt(sum(Psit_full(dirty_image).^2, 2))); 
    mu = nuclear_norm/l21_norm;
    clear dirty_image l21_norm nuclear_norm
    
    % compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
    % SARA space) involved in the reweighting scheme
    B = zeros(size(X0));
    id_center = floor([Ny, Nx]/2) + 1;
    dirac = zeros(Ny, Nx);
    dirac(id_center) = 1;
    AD = A(dirac);
    be = zeros(nChannels, 1);
    for l = 1:nChannels
        temp = zeros(oy*Ny*ox*Nx, 1);
        z = zeros(oy*Ny*ox*Nx, 1);
        for b = 1:numel(y{l})
            noise = (randn(size(y{l}{b})) + 1i*randn(size(y{l}{b})))/sqrt(2);
            temp(W{l}{b}) = temp(W{l}{b}) + G{l}{b}' * (sqrt(aW{l}{b}) .* noise);
            z(W{l}{b}) = z(W{l}{b}) + G{l}{b}' * (sqrt(aW{l}{b}) .* (sqrt(aW{l}{b}).*(G{l}{b} * AD(W{l}{b}))) );  
        end
        B(:,l) = reshape(At(temp),[Ny*Nx,1]);
        be(l) = max(reshape(At(z),[Ny*Nx,1]));     
    end
    B = B/max(be);
    [~ ,S0,~] = svd(B,'econ');
    sig = std(diag(S0));
    sig_bar = std(sqrt(sum(Psit_full(reshape(B, [Ny, Nx, nChannels])).^2,2)));
    clear B S0 be z temp AD dirac id_center Psitw Psit_full Psi1 Psit1
    
    %% HSI parameter structure sent to the  HSI algorithm
    param_HSI.verbose = 2; % print log or not
    param_HSI.nu0 = 1; % bound on the norm of the Identity operator
    param_HSI.nu1 = 1; % bound on the norm of the operator Psi
    param_HSI.gamma0 = 1;
    param_HSI.gamma = mu;  %convergence parameter L1 (soft th parameter)
    param_HSI.rel_obj = 1e-6; % stopping criterion
    param_HSI.max_iter = 10000; % max number of iterations
    
    param_HSI.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
    param_HSI.adapt_eps_start = 250; % minimum num of iter before stating adjustment
    param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
    param_HSI.adapt_eps_tol_out = 1.01; % tolerance outside the l2 ball
    param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
    param_HSI.adapt_eps_rel_var = 1e-4; % bound on the relative change of the solution
    param_HSI.adapt_eps_change_percentage = (sqrt(5)-1)/2; % the weight of the update w.r.t the l2 norm of the residual data
    
    param_HSI.reweight_alpha = 1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
    param_HSI.reweight_alpha_ff = 0.5;
    param_HSI.total_reweights = -1; %30; % -1 if you don't want reweighting
    param_HSI.sig = sig; % estimate of the noise level in SARA space
    param_HSI.sig_bar = sig_bar; % estimate of the noise level in "SVD" space
    param_HSI.use_reweight_steps = 1; % use reweighting steps
    param_HSI.reweight_rel_var = 1e-6; % criterion for performing reweighting
    
    param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid   
    
    %% faceted HyperSARA
    param_HSI.nu2 = Anorm; % upper bound on the norm of the measurement operator A*G
    param_HSI.num_workers = Qx*Qy + ncores_data;

    disp('Faceted HyperSARA')
    disp('-----------------------------------------')

    % spectral tesselation (non-overlapping)
    rg_c = domain_decomposition(ncores_data, channels(end));
    cell_c_chunks = cell(ncores_data, 1);
    y_spmd = cell(ncores_data, 1);
    epsilon_spmd = cell(ncores_data, 1);
    aW_spmd = cell(ncores_data, 1);
    W_spmd = cell(ncores_data, 1);
    G_spmd = cell(ncores_data, 1);

    for i = 1:ncores_data
        cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
        y_spmd{i} = y(cell_c_chunks{i});
        epsilon_spmd{i} = epsilons(cell_c_chunks{i});
        aW_spmd{i} = aW(cell_c_chunks{i});
        W_spmd{i} = W(cell_c_chunks{i});
        G_spmd{i} = G(cell_c_chunks{i});
    end

    param_HSI.ind = 0;

    clear y epsilon aW W G
    %%

    switch algo_version            
        case 'standard' 
            % reference version, overlap between the facets underlying the nuclear norms
            [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                facetHyperSARA(y_spmd, epsilon_spmd, ...
                A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, ncores_data, ...
                wlt_basis, L, nlevel, cell_c_chunks, channels(end), 'whatever.mat', image_name); % [10/10/2019] ok

        case 'weighted' 
            % (piece-wise constant weights, following Audrey's suggestion)
            [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                facetHyperSARA_weighted(y_spmd, epsilon_spmd, ...
                A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, ncores_data, ...
                wlt_basis, L, nlevel, cell_c_chunks, channels(end), 'whatever.mat', image_name); % [10/10/2019] ok

        case 'cst' 
            % cst overlap for the facets undelrying the nuclear norms (d < overlap from sdwt2)
            % d = (power(2, nlevel)-1)*(max(L(:))-2)/2; % use of db8
            [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                facetHyperSARA_cst_overlap(y_spmd, epsilon_spmd, ...
                A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, ncores_data, ...
                wlt_basis, L, nlevel, cell_c_chunks, channels(end), overlap_size, 'whatever.mat', image_name); % [10/10/2019] ok

        case 'cst_weighted' 
            % same as spmd_sct, weight correction (apodization window in this case)
            [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                facetHyperSARA_cst_overlap_weighted(y_spmd, epsilon_spmd, ...
                A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, ncores_data, ...
                wlt_basis, L, nlevel, cell_c_chunks, channels(end), overlap_size, window_type, 'whatever.mat', image_name); % [10/10/2019] ok

            %[xsol,param_HSI,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter] = ...
            %    facetHyperSARA_cst_overlap_weighted_real_data(y_spmd, epsilon_spmd, ...
            %    A, At, aW_spmd, G_spmd, W_spmd, param_HSI, Qx, Qy, Qc2, wlt_basis, L, ...
            %    nlevel, cell_c_chunks, channels(end), d, window_type, 'whatever.mat'); % [10/10/2019] ok basic debugging

        case 'standard2' 
            % alternative implementation (gather primal variables and data on the same nodes)
            % gather image and data on the same nodes (extra communications compared to spmd4 for reweigthing and monitoring variables)
            [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                facetHyperSARA2(y_spmd, epsilon_spmd, ...
                A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, ncores_data, ...
                wlt_basis, L, nlevel, cell_c_chunks, channels(end), 'whatever.mat', image_name); % [10/10/2019] ok

        case 'no_overlap' 
            % same as spmd4, but no overlap for the facets on which the nuclear norms are taken
            [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                facetHyperSARA_no_overlap(y_spmd, epsilon_spmd, ...
                A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, ncores_data, ...
                wlt_basis, L, nlevel, cell_c_chunks, channels(end), 'whatever.mat', image_name); % [10/10/2019] ok

        otherwise
            error('Unknown solver version.')
    end

    time_iter_average = mean(end_iter);
    disp(['snr_x: ', num2str(snr_x)]);
    disp(['asnr_x: ', num2str(snr_x_average)]);
    disp(['Average time per iteration: ', time_iter_average]);

    mkdir('results/')
    save(fullfile(results_path, results_name),'-v7.3','xsol', 'X0', ...
        'param', 'epsilon', 'rel_fval', 'nuclear', 'l21', 'norm_res_out', ...
        'end_iter', 'time_iter_average', 'snr_x', 'snr_x_average');
    fitswrite(xsol,fullfile(results_path, strcat('x_fhs_', algo_version, ...
        '_', window_type, ...
        '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc), ...
        '_overlap=', num2str(overlap_size), ...
        '.fits')))
    
end
