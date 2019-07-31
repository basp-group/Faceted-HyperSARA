clear all
close all
clc

rw = 0;
gamma = 1e-5;
job_id = 0;

flag_algo = 2; tot = 20; Qx = 2; Qy = 1; Qc2 = 1;
parallel_version = '';

addpath fouRed/
addpath lib/
addpath lib/generate_data/
addpath lib/operators/
addpath lib/nufft/
addpath lib/utils/
addpath lib/CubeHelix/
addpath lib/Proximity_operators/code/matlab/indicator/
addpath lib/Proximity_operators/code/matlab/multi/
addpath sdwt2/
addpath src/
addpath src/spmd
addpath src/spmd/dr/
addpath data/
addpath real_data_dr

extract_raw_data = 1;
reduce_data = 1;
generate_eps_nnls = 0;
solve_minimization = 1;
solve_HS = 1;
compute_Anorm = 1;

% select algorithm parameters
window_type = 'triangular';
bool_weights = true; % for the spmd4_new version (50% overlap version)
d = 256; %(power(2, nlevel)-1)*(2*8 - 2); % amount of overlap (absolute number of pixels in the overlap)

% reduction parameters
usingReduction = 1;
usingReductionPar = 0;
usingPrecondition = 1;
normalize_data = 0;
klargestpercent = 20;

%% Config parameters
if compute_Anorm || ~extract_real_data
    Nx = 2560;
    Ny = 1536;
    N = Nx * Ny;
    ox = 2; % oversampling factors for nufft
    oy = 2; % oversampling factors for nufft
    Kx = 8; % number of neighbours for nufft
    Ky = 8; % number of neighbours for nufft
end

%% Extract / load real data (contains Fourier reduction and NNLS)
extract_real_data_dr

%% Compute MAP estimator
if solve_minimization     
    Solver_fouRed_real_data
end
