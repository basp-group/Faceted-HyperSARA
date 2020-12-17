clc; clear all; close all;

addpath lib/generate_data/
addpath lib/operators/
addpath lib/nufft/
addpath lib/utils/
addpath lib/sdwt2/
addpath data/
addpath src/

% Compute residual image (if missing from a given run)
Nx = 1024;
Ny = Nx;
nchannels = 20;
f = linspace(1, 2, nchannels);

coverage_path = 'data/uv_model_N=1024_p=05.fits';
results_path = 'results/W28_1024_m39';
load(fullfile(results_path,'fhs_cw_triangular_N=1024_L=20_p=0.5_Qx=4_Qy=4_Qc=1_overlap=256_snr=60_0_5.6445e-05_1.mat'), 'xsol', 'param');
load('data/y_N=1024_L=20_p=0.5_snr=60.mat', 'y')
y_ = y([1, nchannels]);
xsol = xsol(:,:,[1, nchannels]);
clear y

%% Load uv coverage
uvw = fitsread(coverage_path);
u = uvw(:, 1);
v = uvw(:, 2);
% u = u*f(1)/f(end);
% v = v*f(1)/f(end);
disp('Coverage loaded successfully')
clear uvw

%% Set measurement operator
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

channels = [1, nchannels];
for i = 1:2
    uw = (f(channels(i))/f(1)) * u;
    vw = (f(channels(i))/f(1)) * v;
    
    % compute uniform weights (sampling density) for the preconditioning
    aWw = util_gen_preconditioning_matrix(uw, vw, param_precond);
    
    % set the weighting matrix, these are the natural weights for real data
    nWw = ones(length(uw), 1);
    
%     % set the blocks structure
%     [u1, v1, ~, ~, aW{i}, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    [A, At, G{i}, W{i}] = op_p_nufft([v1 u1], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
end

%% Free memory
clear u v u1 v1 uw vw aWw nW nWw param_block_structure param_precond;

%% Compute residual image
res = zeros(Ny, Nx, 2);
for l = 1:2
    z = zeros(oy*Ny*ox*Nx, 1);
    Fx = A(xsol(:,:,l));
    for b = 1:param_block_structure.equal_partitioning_no
        z(W{l}{b}) = z(W{l}{b}) + G{l}{b}' * (y_{l}{b} - (G{l}{b} * Fx(W{l}{b})));  
    end
    res(:,:,l) = At(z); 
end

%% Save residual image
save('results/W28_1024_m39/res_fhs_cw_triangular_N=1024_L=20_p=0.5_Qx=4_Qy=4_Qc=1_overlap=256_snr=60_0_5.6445e-05_1.mat', '-v7.3', 'res')