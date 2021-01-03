function get_ground_truth(path_uv)

% debugging
% path_uv = '../data';

% clc; clear all; close all;
addpath ../lib/generate_data
addpath ../lib/operators
addpath ../lib/utils
addpath ../lib/measurement-operator/nufft/
seed = 1;

%% Default config parameters
% gridding parameters
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

% preconditioning parameters
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
param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi/4), 0, icdf('norm', 0.75, 0, pi/4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_block_structure.use_manual_partitioning = 0;

%% Re-generate ground truth image cube for spectral faceting experiments
% frequency bandwidth from 1 to 2 GHz
nChannels = 100;
reference_cube_path = '../data/W28_256.fits';
f = linspace(1,2,nChannels);
emission_lines = 0; % insert emission lines on top of the continuous spectra
[x0,~] = Generate_cube_W28(reference_cube_path,f,emission_lines);
[Ny, Nx, ~] = size(x0);
filename = fullfile(path_uv,'vla_7.95h_dt10s.uvw256.mat');
rng(seed);

% recompute operator norm for each channel independently (i.e., values used for SARA)
% load uv-coverage...
param_precond.N = Nx * Ny; % number of pixels in the image
param_precond.Nox = ox*Nx; % number of Fourier points (oversampled plane)
param_precond.Noy = oy*Ny;
load(filename)

u1 = uvw(:, 1);
v1 = uvw(:, 2);
r = sqrt(u1.^2 + v1.^2);
size(r(r>pi));
bmax = max(r);
v1 = v1 * pi/(bmax * 1);
u1 = u1 * pi/(bmax * 1);
u = u1/2;
v = v1/2; 
clear u1 v1 uvw r

operatorNorm = zeros(nChannels, 1);
for l = 1:nChannels
    uw = (f(l)/f(1)) * u;
    vw = (f(l)/f(1)) * v;
    
    % compute uniform weights (sampling density) for the preconditioning
    aWw = util_gen_preconditioning_matrix(uw, vw, param_precond);
    
    % set the weighting matrix, these are the natural weights for real data
    nWw = ones(length(uw), 1);
    
    % set the blocks structure
    [u1, v1, ~, ~, aW, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    [A, At, G, W] = op_p_nufft([v1 u1], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
    F = afclean( @(x) HS_forward_operator_G(x, {G}, {W}, A));
    Ft = afclean( @(y) HS_adjoint_operator_G(y, {G}, {W}, At, Ny, Nx));
    operatorNorm(l) = pow_method_op(F, Ft, [Ny Nx 1]);
end

save('ground_truth_spectral_faceting.mat', '-v7.3', 'x0', 'f', 'operatorNorm');

%% Re-generate ground truth image cube for spectral faceting experiments
% frequency bandwidth from 1 to 2 GHz
nChannels = 20;
reference_cube_path = '../data/W28_1024.fits';
f = linspace(1,2,nChannels);
emission_lines = 0; % insert emission lines on top of the continuous spectra
[x0,~] = Generate_cube_W28(reference_cube_path,f,emission_lines);
[Ny, Nx, ~] = size(x0);
filename = fullfile(path_uv,'vla_7.95h_dt10s.uvw.mat');
rng(seed);

% recompute operator norm for each channel independently (i.e., values used for SARA)
param_precond.N = Nx * Ny; % number of pixels in the image
param_precond.Nox = ox*Nx; % number of Fourier points (oversampled plane)
param_precond.Noy = oy*Ny;
load(filename)

u1 = uvw(:, 1);
v1 = uvw(:, 2);
r = sqrt(u1.^2 + v1.^2);
size(r(r>pi));
bmax = max(r);
v1 = v1 * pi/(bmax * 1);
u1 = u1 * pi/(bmax * 1);
u = u1/2;
v = v1/2; 
clear u1 u2 uvw r

operatorNorm = zeros(nChannels, 1);
for l = 1:nChannels
    uw = (f(l)/f(1)) * u;
    vw = (f(l)/f(1)) * v;
    
    % compute uniform weights (sampling density) for the preconditioning
    aWw = util_gen_preconditioning_matrix(uw, vw, param_precond);
    
    % set the weighting matrix, these are the natural weights for real data
    nWw = ones(length(uw), 1);
    
    % set the blocks structure
    [u1, v1, ~, ~, aW, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    [A, At, G, W] = op_p_nufft([v1 u1], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
    F = afclean( @(x) HS_forward_operator_G(x, {G}, {W}, A));
    Ft = afclean( @(y) HS_adjoint_operator_G(y, {G}, {W}, At, Ny, Nx));
    operatorNorm(l) = pow_method_op(F, Ft, [Ny Nx 1]);
end

save('ground_truth_spatial_faceting.mat', '-v7.3', 'x0', 'f', 'operatorNorm');
