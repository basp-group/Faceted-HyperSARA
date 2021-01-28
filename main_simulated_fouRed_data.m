function main_simulated_fouRed_data(file_name, nChannels, Qx, Qy, Qc, Qc2, ind,...
    flag_algo, window_type, parallel_version, ...
    usingBlocking, usingPrecondition, normalize_data, reduction_version, enable_klargestpercent, fouRed_gamma,...
    flagGenerateGroundtruth, flagUndersampledCube, flagGenerateUv)
% Main script to run the faceted HyperSARA approach with dimensionality reduction on synthetic data
%

addpath fouRed/
addpath lib/
addpath lib/generate_data/
addpath lib/operators/
addpath lib/nufft/
addpath lib/utils/
addpath lib/CubeHelix/
addpath sdwt2/
addpath src/
addpath src/spmd/
addpath src/spmd/dr/
addpath src/spmd/weighted
addpath src/spmd/standard/
addpath data/

% flag_algo = 2; tot = 20; Qx = 2; Qy = 2; Qc = 1; ind = 1; img_size = 512; %2048;
% 
% Qc2 = 5;
% hrs = 6;
% parallel_version = 'spmd4';

% select algorithm parameters
% window_type = 'hamming';
% parallel_version = 'spmd4';
% bool_weights = true; % for the spmd4_new version (50% overlap version)

if strcmp(parallel_version, 'spmd4_cst_weighted') || strcmp(parallel_version, 'spmd4_cst')
    nlevel = 4;
    d = (power(2, nlevel)-1)*(2*8 - 2); % assuming db8 largest wavelet filter
end

% percentage = 0.3;
% 
% usingReduction = 1;
% usingReductionPar = 0;
% usingBlocking = 1;
% usingPrecondition = 1;
% normalize_data = 0;
% reduction_version = 2;      % 1: F\Phi^t 2: G^t (recommend)
% enable_klargestpercent = 0;
% fouRed_gamma = 10;

% Fourier reduction parameters
param_fouRed.enable_klargestpercent = enable_klargestpercent;
% param_fouRed.klargestpercent = 20;
param_fouRed.enable_estimatethreshold = ~enable_klargestpercent;
param_fouRed.gamma = fouRed_gamma;  % reduction parameter


if param_fouRed.enable_estimatethreshold    
    fprintf('Threshold level: remove %d percentile\n', param_fouRed.gamma);
    prob = 1 - param_fouRed.gamma/100;
end

simple_test_1B = 0;

generate_data = 1;

gaussian_UV = 1;
realistic_UV = 0;

flagGenerateUv = 1;
flagSaveUv = 0;

save_data = 0;
load_data = 0;

save_fouRed = 1;
load_fouRed = 0;

save_full_operator = 0;
load_full_operator = 0;

free_memory = 1;
free_Gw = 1;

flagGenerateMeasurements = 1;
flagSaveMeasurements = 0;

solve_minimization = 1;

% % choose the solver you want to use
solve_HS = 1; % wide-band solver (rwLRJS) applied on all the channels
solve_minimization_UQ_HS = 0;

solve_1B = 0; % single-channel solver (SARA) applied for each channel separately
solve_minimization_UQ_1B = 0;

% load_xmap = 1;
% save_xmap = 0;
% 
% load_mask = 1;
% save_mask = 0;

compute_Anorm = 1;
load_Anorm_preconditioned = 0;
save_Anorm_preconditioned = 0;
compute_Anorm_ch = 0;
generate_eps_nnls = 0;

load_Anorm = 0;
save_Anorm = 0;

load_temp_result = 0;

%-- SEE WHERE CHUNK IS USED
% param_HSI.num_workers = 24;
% 
% chunk_width = 2;
% step = 1;
% 
% chunks = [];
% for  i = 1 :step: 15-chunk_width+1
%     chunks = [chunks; i i+chunk_width-1];
% end
% 
% chunks

%% Generate or Load ground-truth
if flagGenerateGroundtruth
%     file_name = ['W28_', num2str(img_size), '.fits'];
%     Nx=img_size;
%     Ny=img_size;
    f = linspace(1,2,nChannels); % frequency bandwidth from 1 to 2 GHz with 15 channels
    emission_lines = 1; % insert emission lines on top of the continuous spectra
    [x0, X0] = Generate_cube(file_name,f,emission_lines);
    [Ny, Nx, ~] = size(x0);
    if flagUndersampledCube
        % new desired dimensions
        % undersampling factor for the channels
        unds = 3; % take 1/unds images
        [x0,X0,f,nChannels] = Generate_undersampled_cube_new(x0,f,Ny,Nx,c,unds);
    end
    fitswrite(x0, cube_path)
    fitsdisp(cube_path)
else
    x0 = fitsread(cube_path);
    [Ny, Nx, ~] = size(x0);
    X0 = reshape(x0, Ny*Nx, nChannels);
end

%% Generate facets
% spectral faceting (interlaved sampling)
id = interleaved_facets(Qc, nChannels);

% create the complete tessellation (spectral)
if ind > 0
    x0 = x0(:,:,id{ind});
    nChannels = size(x0,3);
    f = f(id{ind});
    X0 = reshape(x0, Nx*Ny, nChannels);
end

param_HSI.num_workers = Qx*Qy*2+1;
% %
% if flag_algo == 0 % L11
%     param_HSI.num_workers = 34;
% elseif flag_algo == 1 % HyperSARA
%     param_HSI.num_workers = 1;
% elseif flag_algo == 2 % Faceted HyperSARA
%     param_HSI.num_workers = Qx*Qy*2+1;  %%%%%%%%%%% TO BE SET BY P.-A.
% end
% param_HSI.num_workers

%% Generate or Load subsampling mask
percentage = 0.3; %[0.3 0.05 0.1 0.02 0.01 0.2]; % sampling rate

%% TO BE MODIFIED FROM HERE

%% Generate or Load subsampling mask
%% config parameters
N = Nx * Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

%% preconditioning parameters
param_precond.N = N; % number of pixels in the image
param_precond.Nox = ox*Nx; % number of pixels in the image
param_precond.Noy = oy*Ny; % number of pixels in the image
param_precond.gen_uniform_weight_matrix = 1; %set weighting type
param_precond.uniform_weight_sub_pixels = 1;


%% block structure

param_block_structure.use_density_partitioning = 0;
param_block_structure.density_partitioning_no = 1;

param_block_structure.use_uniform_partitioning = 0;
param_block_structure.uniform_partitioning_no = 4;

param_block_structure.use_equal_partitioning = 1;
param_block_structure.equal_partitioning_no = 1;

param_block_structure.use_manual_frequency_partitioning = 0;
param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi/4), 0, icdf('norm', 0.75, 0, pi/4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_block_structure.use_manual_partitioning = 0;

%% samplig pattern parameters
% options 'gaussian', 'file', 'gaussian+large-holes', 'file+undersample', ''gaussian+missing-pixels'
sampling_pattern = 'gaussian+large-holes';
sparam.p = percentage; % number of measurements as proportion of number of pixels to recover
sparam.hole_number = 8000; % number of holes to introduce for 'gaussian+large-holes'
sparam.hole_prob = 0.05; % probability of single pixel hole for 'gaussian+missing-pixels'
sparam.hole_size = pi/60; % size of the missing frequency data
sparam.fpartition = [pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [0, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.25*pi, 0, 0.25*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
sparam.sigma = pi/3; % variance of the gaussion over continous frequency
sparam.sigma_holes = pi/3; % variance of the gaussion for the holes

%% generate the sampling pattern
sparam.N = N; % number of pixels in the image
sparam.Nox = ox*Nx; % number of pixels in the image
sparam.Noy = oy*Ny; % number of pixels in the image


%%
if flagGenerateUv
    [u, v, ~, ~, ~] = util_gen_sampling_pattern(sampling_pattern, sparam);
    % remove all the visibilities outside the range [-pi, pi]
    r = sqrt(u{1}.^2 + v{1}.^2);

    u1 = u{1};
    v1 = v{1};
    [mm] = find(r <= pi);
    
    u1 = u1(mm);
    v1 = v1(mm);

    r = sqrt(u1.^2 + v1.^2);
    size(r(r>pi))

    u1 = u1/2;
    v1 = v1/2;
    
    if flagSaveUv
        save(['./data/uv_' num2str(c) '_HI_final_' num2str(percentage(k)) '.mat'],'-v7.3', 'u', 'v');
    end
else
    load(['./data/uv_' num2str(c) '_HI_final_' num2str(percentage(k)) '.mat']);
    disp('coverage loaded successfully ;)')
end

%scatter(u1,v1,'r.')

%%
uw = cell(nChannels, 1);
vw = cell(nChannels, 1);
aW = cell(nChannels, 1);
G = cell(nChannels, 1);
W = cell(nChannels, 1);
for i = 1:nChannels
    
    uw{i} = (f(i)/f(1)) * u1;
    vw{i} = (f(i)/f(1)) * v1;

    %% compute uniform weights (sampling density) for the preconditioning
    [aWw] = util_gen_preconditioning_matrix(uw{i}, vw{i}, param_precond);
    
    % set the weighting matrix, these are the natural weights for real data
    nWw = ones(length(uw{i}), 1);
    
    % set the blocks structure
    [u, v, ~, ~, aW{i}, nW] = util_gen_block_structure(uw{i}, vw{i}, aWw, nWw, param_block_structure);
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    
    [A, At, G{i}, W{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
    
end

%% Free memory
clear u v u1 v1 uw vw aWw nW nWw r mm uvidx

%%
sigma_noise = 0.1;
input_snr = 40; % noise level (on the measurements)
% Generate or Load subsampling mask

if flagGenerateMeasurements
    [y, G, ~, ~, ~] = Generate_fouRed_Measurements(x0, G, A, W, input_snr, sigma_noise, usingBlocking, normalize_data);
    if flagSaveMeasurements
        save(['./data/data_' num2str(c) '_HI_final_' num2str(percentage(k)) '.mat'],'-v7.3', 'y', 'G', 'W', 'noise', 'epsilon', 'epsilons');
    end
else 
    load(['./data/data_' num2str(c) '_HI_final_' num2str(percentage(k)) '.mat'], 'y', 'G', 'W', 'noise', 'epsilon', 'epsilons')
end
    
%% Fourier reduction

[yT, H, W, T, aW, Wm, epsilon, epsilons] = Fourier_reduction_ch(y, G, reduction_version, fouRed_gamma, param_fouRed, save_fouRed, load_fouRed);

clear y G
%% Faceted HyperSARA with DR

Solver_simulated_fouRed_data_FINAL_clean
