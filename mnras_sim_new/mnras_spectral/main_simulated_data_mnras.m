function main_simulated_data_mnras(image_name, nChannels, Qx, Qy, Qc, ...
    algo_version, window_type, ncores_data, ind, overlap_fraction, nReweights, ...
    flag_generateCube, flag_generateVisibilities, ...
    flag_computeOperatorNorm, flag_solveMinimization, ...
    cube_path, coverage_path, gam, rw, flag_homotopy, ... 
    flag_computeLowerBounds, rwtype, gam_bar, exp_type, superresolution_factor)
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
%        - 'cw'            cst_weighted: constant overlap taken for the 
%                          faceted nuclear norm, using spatial weights 
%                          (apodization window)
%        - 'no'            no_overlap: no overlap for the faceted nuclear 
%                          norm prior
%     window_type (string): type of apodization window considered for the 
%                           faceted nuclear norm prior. Only active with 
%                           the following versions of the algorithm:
%        - 'triangular'
%        - 'hamming'
%        - 'pc' (piecewise-constant)
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

%% PARAMETERS FOR DEBUGGING

% image_name = 'W28_512'; %'cygASband_Cube_H'; %'W28_512';
% exp_type = 'local_test'; % 'spectral', 'spatial', 'test'

% Qx = 2; % 4
% Qy = 1; % 4
% Qc = 1;
% nReweights = 1;
% algo_version = 'sara'; % 'cw', 'hypersara', 'sara';
% window_type = 'triangular'; % 'hamming', 'pc'
% flag_generateVisibilities = 0;
% flag_computeOperatorNorm = 0;
% flag_computeLowerBounds = 0;
% flag_solveMinimization = true;
% ncores_data = 1; % number of cores assigned to the data fidelity terms (groups of channels)
% ind = 1; % index of the spectral facet to be reconstructed
% gam = 1;
% gam_bar = 1;
% coverage_path = "data/msSpecs.mat"; % "data/vla_7.95h_dt10s.uvw256.mat";

% rw = 1;
% rwtype = 'dirty'; % ground_truth, heuristic
% flag_homotopy = 1;
% overlap_fraction = 0.5;

% nChannels = 5;
% flag_generateCube = 1;
% cubepath = @(nchannels) strcat(image_name, '_L', num2str(nchannels));
% cube_path = cubepath(nChannels);
% flag_generateCoverage = 0;
% flag_generateUndersampledCube = 0; % Default 15 channels cube with line emissions

%%

% fixed parameters (in the mnras experiments)
flag_generateUndersampledCube = false;
flag_generateCoverage = false;
p = 1;
seed = 1;
rng(seed);
% sigma_noise = 0.1
% T = 1500; % to be set depending on the value of p
% hrs = 6;
% kernel = 'minmax:tuned'; % 'kaiser', 'minmax:tuned'
generate_eps_nnls = false;
save_data = true;
speed_of_light = 299792458;

%%
format compact;

disp('MNRAS configuration')
disp(['Algorithm version: ', algo_version]);
disp(['Reference image: ', image_name]);
disp(['nchannels: ', num2str(nChannels)]);
disp(['Number of facets Qy x Qx : ', num2str(Qy), ' x ', num2str(Qx)]);
disp(['Number of spectral facets Qc : ', num2str(Qc)]);
disp(['Overlap fraction: ', strjoin(strsplit(num2str(overlap_fraction)), '_')]);
% disp(['Number of data points p per frequency (as a fraction of image size): ', num2str(p)]);
% disp(['Input SNR: ', num2str(input_snr)]);
disp(['Generating image cube: ', num2str(flag_generateCube)]);
disp(['Generating coverage: ', num2str(flag_generateCoverage)]);
disp(['Generating visibilities: ', num2str(flag_generateVisibilities)]);

addpath ../../lib/generate_data/
addpath ../../lib/operators/
addpath ../../lib/measurement-operator/nufft/
addpath ../../lib/measurement-operator/lib/operators/
addpath ../../lib/utils/
addpath ../../lib/faceted-wavelet-transform/src
addpath ../../data/
addpath ../../src_mnras/
if strcmp(algo_version, "hypersara")
    addpath ../../src_mnras/serial
else
    addpath ../../src_mnras/spmd
    addpath ../../src_mnras/spmd/weighted
end

% setting paths to results and reference image cube
% coverage_path = strcat(coverage_path, '.fits');
data_path = '../../data/';
results_path = fullfile('results/', strcat(image_name, '_', exp_type));
mkdir(data_path)
mkdir(results_path)

%% Generate/load ground-truth image cube
% reference_cube_path = strcat(data_path, image_name, '.fits');
% if flag_generateCube
%     % frequency bandwidth from 1 to 2 GHz
%     f = linspace(1,2,nChannels);
%     emission_lines = 0; % insert emission lines on top of the continuous spectra
%     % [x0,X0] = Generate_cube(reference_cube_path,f,emission_lines);
%     [x0,X0] = Generate_cube_W28(reference_cube_path,f,emission_lines);
%     [Ny, Nx, nChannels] = size(x0);
%     if flag_generateUndersampledCube
%         % undersampling factor for the channels
%         unds = 4; % take 1/unds images
%         [x0,X0,f,nChannels] = Generate_undersampled_cube(x0,f,Ny,Nx,nChannels,unds);
%     end
%     fitswrite(X0.', strcat(cube_path, '.fits'));
%     fitsdisp(strcat(cube_path, '.fits'));
% else
%     X0 = fitsread(strcat(cube_path, '.fits')).';
%     Nx = sqrt(size(X0, 1));
%     Ny = Nx;
%     nChannels = size(X0, 2);
%     x0 = reshape(X0, [Ny, Nx, nChannels]);
% end
% % frequency bandwidth from 1 to 2 GHz
% f = linspace(1, 2, nChannels);

%%
%! see how to load the cube prepared by Arwa
%! load the appropriate portion of the reference image cube
% if spectral faceting, just load the intersting portion of the full image cube
switch exp_type
    case "spatial"
        image_name = 'cygASband_Cube_H';
        spectral_downsampling = 5;
        spatial_downsampling = 1;     
    case "spectral"
        image_name = 'cygASband_Cube_L';
        spectral_downsampling = 1;
        spatial_downsampling = 1;
    case "test"
        image_name = 'cygASband_Cube_1024x512x20';
        spectral_downsampling = 1;
        spatial_downsampling = 1;
    case "local_test"
        image_name = 'cygASband_Cube_1024x512x20';
        spectral_downsampling = 5;
        spatial_downsampling = 2;
        coverage_path = "data/vla_7.95h_dt10s.uvw256.mat";
    case "old_local_test"
        image_name = 'cubeW28';
        spectral_downsampling = 20;
        spatial_downsampling = 4;
        coverage_path = "data/vla_7.95h_dt10s.uvw256.mat";
    otherwise
        error("Unknown experiment type")
end

reference_cube_path = fullfile(data_path, strcat(image_name, '.fits'));
info        = fitsinfo(reference_cube_path);
rowend      = info.PrimaryData.Size(1);
colend      = info.PrimaryData.Size(2);
sliceend    = info.PrimaryData.Size(3);
if strcmp(algo_version, 'sara')
    x0 = fitsread(reference_cube_path, 'primary',...
              'Info', info,...
              'PixelRegion',{[1 spatial_downsampling rowend], ...
              [1 spatial_downsampling colend], ...
              ind});
else
    x0 = fitsread(reference_cube_path, 'primary',...
              'Info', info,...
              'PixelRegion',{[1 spatial_downsampling rowend], ...
              [1 spatial_downsampling colend], ...
              [1 spectral_downsampling sliceend]});
end
nChannels = floor(sliceend/spectral_downsampling);
clear info rowend colend

[Ny, Nx, nchans] = size(x0);
N = Nx*Ny;
X0 = reshape(x0, [N, nchans]);
input_snr = 50*ones(nchans, 1); % input SNR (in dB)

% frequency used to generate the 2 reference cubes
nu0 = 2.0525e9; % starting freq
dnu = 16e6;    % freq step
% superresolution_factor = 2; %? keep as is, or change order of magnitude
% pixelsizeH = 0.2*1024/3000;
% pixelsizeL = 0.2*1024/512;
f = nu0 + (0:spectral_downsampling:(sliceend-1))*dnu;

%% Generate spectral facets (interleaved sampling)
if strcmp(algo_version, 'sara')
    window_type = 'none';
    Qc = nChannels; %! handle each channel separately
    Qx = 1;
    Qy = 1;
elseif strcmp(algo_version, 'hypersara')
    window_type = 'none';
    Qx = 1;
    Qy = 1;
end

overlap_size = get_overlap_size([Ny, Nx], [Qy, Qx], overlap_fraction);

id = split_range_interleaved(Qc, nChannels);
if Qc > 1 && ind > 0 && ~strcmp(algo_version, 'sara')
    x0 = x0(:,:,id{ind});
    nchans = size(x0,3);
    f = f(id{ind}); %? do not do this, otherwise the uv-normalization is
    % completely wrong?
    X0 = reshape(x0,Nx*Ny,nchans);
    input_snr = input_snr(id{ind});
end
channels = 1:nchans;

%% Setup name of results file
data_name_function = @(nchannels) strcat('y_', ...
    exp_type,'_',image_name, '_srf=', num2str(superresolution_factor), ...
    '_Ny=',num2str(Ny),'_Nx=',num2str(Nx),'_L=', num2str(nchannels),'.mat');

results_name_function = @(nchannels) strcat(exp_type,'_',image_name,'_', ...
    algo_version,'_',window_type, '_srf=', num2str(superresolution_factor), ...
    '_Ny=',num2str(Ny), '_Nx=',num2str(Nx), '_L=',num2str(nchannels), ...
    '_Qy=', num2str(Qy), '_Qx=', num2str(Qx), '_Qc=', num2str(Qc), ...
    '_ind=', num2str(ind), '_gam=', num2str(gam), '_gambar=', num2str(gam_bar), ...
    '_overlap=', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
    '_rw_type=', rwtype, ...
    '.mat');

temp_results_name = @(nchannels) strcat(exp_type,'_',image_name,'_', ...
    algo_version,'_',window_type, '_srf=', num2str(superresolution_factor), ...
    '_Ny=',num2str(Ny), '_Nx=',num2str(Nx), '_L=',num2str(nchannels), ...
    '_Qy=', num2str(Qy), '_Qx=', num2str(Qx), '_Qc=', num2str(Qc), ...
    '_ind=', num2str(ind), '_gam=', num2str(gam), '_gambar=', num2str(gam_bar), ...
    '_overlap=', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
    '_rw_type=', rwtype);

warm_start = @(nchannels) strcat(temp_results_name(nchannels),'_rw=', num2str(rw), '.mat');

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
    coverage_path	
    load(coverage_path, 'uvw');
    size(uvw)
    
    %! Abdullah's version: normalization w.r.t. the minimum frequency
    % u1 = uvw(:, 1);
    % v1 = uvw(:, 2);
    % r = sqrt(u1.^2 + v1.^2);
    % bmax = max(r);
    % v1 = v1 * pi/(bmax * 1); 
    % u1 = u1 * pi/(bmax * 1);
    % u = u1/2;
    % v = v1/2;
    
    %! normalize u,v coverage w.r.t. the highest frequency (i.e., uv expressed in
    % units of the smallest wavelenght, associated with the highest frequency)
    u1 = uvw(:, 1)*f(end)/speed_of_light;
    v1 = uvw(:, 2)*f(end)/speed_of_light;
    bmax = max(sqrt(u1.^2 + v1.^2));
    % cellsize = 3600*180/(superresolution_factor*2*pi*bmax); % in arcsec
    u = u1*pi/(superresolution_factor*2*bmax);
    v = v1*pi/(superresolution_factor*2*bmax);
    size(u)
    disp('Coverage loaded successfully')
    clear uvw u1 v1
end

% setup measurement operator
for i = 1:nchans
    %! previous version
    % uw = (f(i)/f(1)) * u;
    % vw = (f(i)/f(1)) * v;

    %? need to fix the frequencies? (just select the right ones)
    %? i.e., need to normalize w.r.t. max over all available frequencies?
    uw = (f(i)/f(end)) * u;
    vw = (f(i)/f(end)) * v;
    
    % compute uniform weights (sampling density) for the preconditioning
    aWw = util_gen_preconditioning_matrix(uw, vw, param_precond);
    
    % set the weighting matrix, these are the natural weights for real data
    nWw = ones(length(uw), 1);
    
    % set the blocks structure
    [u1, v1, ~, ~, aW{i}, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    [A, At, G{i}, W{i}] = op_p_nufft([v1 u1], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
%     for b = 1:numel(G{i})
%         G{i}{b} = G{i}{b}/sigma_noise; %! add variance normalisation to use the same procedure as for real data to estimate the reweighting lower bounds
%     end
end

%% Free memory
clear u v u1 v1 uw vw aWw nW nWw param_block_structure param_precond;

%% Generate/load visibilities (generate only full spectral dataset)
if flag_generateVisibilities
    param_l2_ball.type = 'sigma';
    param_l2_ball.sigma_ball = 2;
    % [y0, y, Nm, sigma_noise] = util_gen_measurements(x0, G, W, A, input_snr);
    
    %! see if it realy needs to be hard-coded this way!
    % [y0, y, Nm] = util_gen_measurements_sigma(x0, G, W, A, 1, seed);
    % [epsilon,epsilons] = util_gen_data_fidelity_bounds(y, Nm, param_l2_ball, 1); %! sigma_noise modified to account for variance normalisation (Theta) % sigma_noise
    
    %! [15/04/2021] agreed with Arwa and Yves 
    [y0, y, Ml, Nm, sigma_noise] = util_gen_measurements_snr(x0, G, W, A, input_snr,seed);
    [epsilon,epsilons] = util_gen_data_fidelity_bounds2(y, Ml, param_l2_ball, sigma_noise);    
    
    if save_data
        save(fullfile(results_path,data_name), '-v7.3', 'y0', 'y', 'epsilon', 'epsilons', 'sigma_noise');
    end
    clear y0 Nm epsilon;
else
    %! if spectral faceting or SARA, only load the portion of the data
    %! needed for the analysis
    load(fullfile(results_path,data_name), 'y', 'epsilons', 'sigma_noise');
    %! need to completely restructure data format to avoid reloading
    %! everything
    %m = matfile(fullfile(data_path,data_name));
    %y = cell(numel(id{ind}),1);
    %epsilons = cell(numel(id{ind}),1);
    %for k = 1:numel(id{ind})
    %    y(k) = cell2mat(m.y(id{ind},1));
    %    epsilons(k) = m.epsilons(id{ind}(k),1);
    %end
    disp('Data loaded successfully')
end
y = y(id{ind});
epsilons = epsilons(id{ind});

%% Compute operator norm
if strcmp(algo_version, 'sara')
    % Compute measurement operator spectral norm for each channel individually
    if flag_computeOperatorNorm
        %Anorm_ch(i) = pow_method_op(@(x) sqrt(cell2mat(aW{i})) .* (Gw{i}*A(x)), @(x) real(At(Gw{i}' * (sqrt(cell2mat(aW{i})) .* x))), [Ny Nx 1]);
        F = afclean( @(x) HS_forward_operator_precond_G(x, G, W, A, aW));
        Ft = afclean( @(y) HS_adjoint_operator_precond_G(y, G, W, At, aW, Ny, Nx));
        Anorm = op_norm(F, Ft, [Ny Nx nchans], 1e-8, 200, 2);
        save(fullfile(results_path, ...
            strcat('Anorm_', ...
            algo_version, ...
            '_Ny=',num2str(Ny),'_Nx=',num2str(Nx), '_L=',num2str(nChannels), ...
            '_Qc=',num2str(Qc),'_ind=',num2str(ind), '_ch=', num2str(ind), '.mat')),'-v7.3', 'Anorm');
    else
        load(fullfile(results_path, ...
            strcat('Anorm_', ...
            algo_version, ...
            '_Ny=',num2str(Ny),'_Nx=',num2str(Nx), '_L=',num2str(nChannels), ...
            '_Qc=',num2str(Qc),'_ind=',num2str(ind), '_ch=', num2str(ind), '.mat')));
    end
else
    if flag_computeOperatorNorm
        % Compute full measurement operator spectral norm
        F = afclean( @(x) HS_forward_operator_precond_G(x, G, W, A, aW));
        Ft = afclean( @(y) HS_adjoint_operator_precond_G(y, G, W, At, aW, Ny, Nx));
        Anorm = op_norm(F, Ft, [Ny Nx nchans], 1e-8, 200, 2);
        save(fullfile(results_path, ...
            strcat('Anorm_', ...
            algo_version, ...
            '_Ny=',num2str(Ny), '_Nx=',num2str(Nx),'_L=', num2str(nChannels), ...
            '_Qc=',num2str(Qc),'_ind=',num2str(ind), '.mat')),'-v7.3', 'Anorm');
       
        % % operator norm per channel "l"
        % Anorm_channel = zeros(nchans, 1);
        % for l = 1:nchans
        %     F = afclean( @(x) HS_forward_operator_precond_G(x, G(l), W(l), A, aW(l)));
        %     Ft = afclean( @(y) HS_adjoint_operator_precond_G(y, G(l), W(l), At, aW(l), Ny, Nx));
        %     Anorm_channel(l) = op_norm(F, Ft, [Ny Nx 1], 1e-8, 200, 2); 
        % end
        % save(fullfile(results_path,strcat('Anorm_channel_N=',num2str(Nx), ...
        %     '_L=',num2str(nChannels),'_Qc=',num2str(Qc),'_ind=',num2str(ind), '.mat')),'-v7.3', 'Anorm_channel');
    else
        load(fullfile(results_path, ...
            strcat('Anorm_', ...
            algo_version, ...
            '_Ny=',num2str(Ny), '_Nx=',num2str(Nx),'_L=', num2str(nChannels), ...
            '_Qc=',num2str(Qc),'_ind=',num2str(ind), '.mat')));
    end
end
    
%% Generate initial epsilons by performing imaging with NNLS on each data block separately
if generate_eps_nnls
    % param_nnls.im = im; % original image, used to compute the SNR
    param_nnls.verbose = 2; % print log or not
    param_nnls.rel_obj = 1e-5; % stopping criterion
    param_nnls.max_iter = 1000; % max number of iterations
    param_nnls.sol_steps = [inf]; % saves images at the given iterations
    param_nnls.beta = 1;
    % solve nnls per block
    for i = 1:nchans
        eps_b{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            % printf('solving for band %i\n\n',i)
            [~,eps_b{i}{j}] = fb_nnls_blocks(y{i}{j}, A, At, G{i}{j}, W{i}{j}, param_nnls);
        end
    end
    mkdir('data/')
    save('data/eps.mat','-v7.3', 'eps_b');
end

%% Solver
if flag_solveMinimization
    % wavelets
    nlevel = 4; % depth of the wavelet decompositions
    %! always specify Dirac basis ('self') in last position if used in the
    %! SARA dictionary 
    wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; 
    filter_length = [2*(1:8)'; 0]; % length of the filters (0 corresponding to the 'self' basis)
    
    %! -- TO BE CHECKED
    % compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
    % SARA space) involved in the reweighting scheme
    if strcmp(algo_version, 'sara')
        dwtmode('zpd')
        [Psi1, Psit1] = op_p_sp_wlt_basis_fhs(wlt_basis, nlevel, Ny, Nx);
        P = length(Psi1);
        for k = 1 : P
            f = '@(x_wave) HS_forward_sparsity(x_wave,Psi1{';
            f = sprintf('%s%i},Ny,Nx);', f,k);
            Psi{k} = eval(f);
            b(k) = size(Psit1{k}(zeros(Ny,Nx,1)),1);
            ft = ['@(x) HS_adjoint_sparsity(x,Psit1{' num2str(k) '},b(' num2str(k) '));'];
            Psit{k} = eval(ft);
        end
        if flag_computeLowerBounds
            %! to be updated (with the different options)
            [sig, max_psf, mu, l21_norm, l21_norm_x0, dirty_image] = ...
                compute_reweighting_lower_bound_sara(y, W, G, A, At, Ny, Nx, ...
                oy, ox, wlt_basis, filter_length, nlevel, sigma_noise, rwtype, x0, Anorm);
            g = gam;
            gam = gam*mu;
            save(fullfile(results_path, ...
                strcat('lower_bounds_', ...
                algo_version, '_srf=', num2str(superresolution_factor), ...
                '_Ny=',num2str(Ny), '_Nx=',num2str(Nx), '_L=',num2str(nChannels),...
                '_ind=', num2str(ind), ...
                '_gam=', num2str(g), ...
                '_rwtype=', rwtype, '.mat')), ...
                'sig', 'max_psf', 'l21_norm', 'l21_norm_x0', 'dirty_image', 'gam');
        else
            load(fullfile(results_path, ...
                strcat('lower_bounds_', ...
                algo_version, '_srf=', num2str(superresolution_factor), ...
                '_Ny=',num2str(Ny), '_Nx=',num2str(Nx), '_L=',num2str(nChannels),...
                '_ind=', num2str(ind), ...
                '_gam=', num2str(gam), ...
                '_rwtype=', rwtype, '.mat')), ...
                'sig', 'max_psf', 'gam');
        end
    else
        fprintf('Normalization factors gamma = %e, gamma_bar = %e \n', gam, gam_bar);
        if flag_computeLowerBounds
            %! use normalization by operator norm
            [sig, sig_bar, mu0, mu, mu_bar, ~, max_psf, l21_norm, nuclear_norm, l21_norm_x0, nuclear_norm_x0, dirty_image] = compute_reweighting_lower_bound(y, W, G, A, At, Ny, Nx, oy, ox, ...
                nchans, wlt_basis, filter_length, nlevel, sigma_noise, rwtype, algo_version, Qx, Qy, overlap_size, window_type, x0, Anorm);
            fprintf('%s: upsilon_bar min = %e, max = %e\n', rwtype, min(sig_bar), max(sig_bar));
            fprintf('%s: mu_bar = %e, mu = %e\n', rwtype, mu_bar, mu);
            fprintf('upsilon = %e \n', sig);
            
            %! recompute the value for gam (ratio between l21 and nuclear norm)
            % g = gam;
            % gam = gam*nuclear_norm/l21_norm;
            g = gam;
            g_bar = gam_bar;
            gam = gam*mu;
            gam_bar = gam_bar*mu_bar;
            save(fullfile(results_path, ...
                strcat('lower_bounds_', ...
                algo_version, ...
                '_srf=', num2str(superresolution_factor), ...
                '_Ny=',num2str(Ny), '_Nx=',num2str(Nx), '_L=',num2str(nChannels), ...
                '_Qy=', num2str(Qy), '_Qx=', num2str(Qx), '_Qc=', num2str(Qc), ...
                '_ind=', num2str(ind), ...
                '_gam=', num2str(g), '_gambar=', num2str(g_bar), ...
                '_rwtype=', rwtype, '.mat')), ...
                    'sig', 'sig_bar', 'max_psf', 'l21_norm', 'nuclear_norm', 'dirty_image', 'gam', 'l21_norm_x0', 'nuclear_norm_x0', 'gam_bar');
        else
            load(fullfile(results_path, ...
                strcat('lower_bounds_', ...
                algo_version, ...
                '_srf=', num2str(superresolution_factor), ...
                '_Ny=',num2str(Ny), '_Nx=',num2str(Nx), '_L=',num2str(nChannels), ...
                '_Qy=', num2str(Qy), '_Qx=', num2str(Qx), '_Qc=', num2str(Qc), ...
                '_ind=', num2str(ind), ...
                '_gam=', num2str(gam), '_gambar=', num2str(gam_bar), ...
                '_rwtype=', rwtype, '.mat')), ...
                'sig', 'sig_bar', 'max_psf', 'l21_norm', 'nuclear_norm', 'gam', 'gam_bar');
        end
        fprintf('Reweighting type (for SVD log prior): %s \n', rwtype);
        fprintf('nuclear_norm_dirty = %e, l21_norm_dirty = %e\n', l21_norm, nuclear_norm);
        fprintf('gam_bar = %e\n', gam_bar);
        fprintf('sig_bar = %e\n', sig_bar);
    end
    clear dirty_image
    fprintf('%s: gam = %e, sig = %e\n\n', rwtype, gam, sig);
    %! --
    
    %% HSI parameter structure sent to the  HSI algorithm
    param_HSI.verbose = 2; % print log or not
    param_HSI.nu0 = 1; % bound on the norm of the Identity operator
    param_HSI.nu1 = 1; % bound on the norm of the operator Psi
    param_HSI.nu2 = Anorm; % upper bound on the norm of the measurement operator
    param_HSI.gamma0 = gam_bar;  % regularization parameter nuclear norm
    param_HSI.gamma = gam; % regularization parameter l21-norm (soft th parameter) %! for SARA, take the value given as an input to the solver
    param_HSI.cube_id = ind;  % id of the cube to be reconstructed (if spectral faceting active)

    % pdfb
    % param_HSI.pdfb_min_iter = 100; % minimum number of iterations
    param_HSI.pdfb_max_iter = 2000; % maximum number of iterations
    param_HSI.pdfb_rel_var = 1e-5; % relative variation tolerance
    param_HSI.pdfb_fidelity_tolerance = 1.01; % tolerance to check data constraints are satisfied %! this value seems quite stringent in practice
    
    % epsilon update scheme
    param_HSI.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
    param_HSI.adapt_eps_start = 200; % minimum num of iter before stating adjustment
    param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
    param_HSI.adapt_eps_tol_out = 1.01; % tolerance outside the l2 ball
    param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
    param_HSI.adapt_eps_rel_var = 1e-5; % bound on the relative change of the solution
    param_HSI.adapt_eps_change_percentage = (sqrt(5)-1)/2; % the weight of the update w.r.t the l2 norm of the residual data
    
    %! -- TO BE CHECKED
    param_HSI.reweighting_rel_var = 1e-5;       % relative variation (reweighting)
    if flag_homotopy
        param_HSI.reweighting_alpha = 10;
        param_HSI.reweighting_min_iter = 5; % minimum number of reweighting iterations, weights updated reweighting_min_iter times
        param_HSI.reweighting_alpha_ff = (1/param_HSI.reweighting_alpha)^(1/(param_HSI.reweighting_min_iter-1)); % reach the floor level after min_iter-2 updates of the weights
        % 0.63 -> otherwise need 10 reweights minimum
    else
        param_HSI.reweighting_min_iter = 1; % minimum number of reweighting iterations
        param_HSI.reweighting_alpha = 1;
        param_HSI.reweighting_alpha_ff = 1;
    end
    %! --
    param_HSI.reweighting_max_iter = max(nReweights, param_HSI.reweighting_min_iter); % maximum number of reweighting iterations reached, weights updated nReweights times
    param_HSI.reweighting_sig = sig; % estimate of the noise level in SARA space
    if ~strcmp(algo_version, 'sara') %! if HyperSARA or faceted HyperSARA
        % estimate of the noise level in "SVD" spaces
        param_HSI.reweighting_sig_bar = sig_bar; 
    end
    param_HSI.max_psf = max_psf;

    % ellipsoid projection parameters (when preconditioning is active)
    param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid   

    %%
    if strcmp(algo_version, 'sara')
        disp('SARA')
        disp('-----------------------------------------')

        [xsol,param,v1,v2,g,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_val,l11,norm_res,res,t_l11,t_master,end_iter] = ...
            sara2(y, epsilons, A, At, aW, G, W, Psi, Psit, param_HSI, fullfile(results_path,warm_start(nChannels)), ...
            fullfile(results_path,temp_results_name(nChannels)), x0, flag_homotopy, ncores_data); %! in this case, ncores_data corresponds to the number of workers for the wavelet transform (9 maximum)
        
        time_iter_average = mean(end_iter);
        disp(['Average time per iteration: ', num2str(time_iter_average)]);

        mkdir('results/')
        save(fullfile(results_path, results_name),'-v7.3','xsol', 'X0', ...
        'param', 'epsilon', 'rel_val', 'l11', 'norm_res', ...
        'end_iter', 'time_iter_average', 't_l11','t_master', 'res');
        fitswrite(xsol,fullfile(results_path, strcat('x_', image_name, '_', algo_version, ...
        '_srf=', num2str(superresolution_factor), ...
        '_', window_type, ...
        '_Qy=', num2str(Qy), '_Qx=', num2str(Qx), '_Qc=', num2str(Qc), ...
        '_ind=', num2str(ind), ...
        '_homotopy=', num2str(flag_homotopy), ...
        'rwtype=', rwtype, ...
        '.fits')))
    else
        %%
        disp('Faceted HyperSARA')
        disp('-----------------------------------------')

        % spectral tesselation (non-overlapping)
        rg_c = split_range(ncores_data, channels(end));
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
        clear y epsilon aW W G
        
        %%
        switch algo_version 
            case 'hypersara'
                % reference HyperSARA version
                [xsol,param,epsilon,t,rel_val,nuclear,l21,norm_res_out,res,end_iter,snr_x,snr_x_average] = ...
                    hyperSARA2(y_spmd, epsilon_spmd, ...
                    A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, ncores_data, ...
                    wlt_basis, nlevel, cell_c_chunks, channels(end), fullfile(results_path,warm_start(nChannels)), fullfile(results_path,temp_results_name(nChannels)), flag_homotopy);

            case 'cw'
                [xsol,param,epsilon,t,rel_val,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                    facetHyperSARA_cw2(y_spmd, epsilon_spmd, ...
                    A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, ncores_data, wlt_basis, ...
                    filter_length, nlevel, cell_c_chunks, channels(end), overlap_size, window_type, fullfile(results_path,warm_start(nChannels)), fullfile(results_path,temp_results_name(nChannels)), flag_homotopy);
            otherwise
                error('Unknown solver version.')
        end
    end
    
    time_iter_average = mean(end_iter);
    disp(['snr_x: ', num2str(snr_x)]);
    disp(['asnr_x: ', num2str(snr_x_average)]);
    disp(['Average time per iteration: ', num2str(time_iter_average)]);

    mkdir('results/')
    save(fullfile(results_path, results_name),'-v7.3','xsol', 'X0', ...
        'param', 'epsilon', 'rel_val', 'nuclear', 'l21', 'norm_res_out', ...
        'end_iter', 'time_iter_average', 'snr_x', 'snr_x_average');
    fitswrite(xsol,fullfile(results_path, strcat('x_', image_name, '_', algo_version, ...
        '_', window_type, ...
        '_srf=', num2str(superresolution_factor), ...
        '_Qy=', num2str(Qy), '_Qx=', num2str(Qx), '_Qc=', num2str(Qc), ...
        '_ind=', num2str(ind), ...
        '_overlap=', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
        '_homotopy=', num2str(flag_homotopy), ...
        'rwtype=', rwtype, ...
        '.fits')))
end
