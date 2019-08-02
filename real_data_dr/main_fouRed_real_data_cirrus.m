function main_fouRed_real_data_cirrus(gamma, rw, job_id)

% gamma = 1e-5;
flag_algo = 2;
Qx = 5; 
Qy = 3; 
Qc2 = 20;
% parallel_version = '';

extract_real_data = false;
reduce_data = 1;
generate_eps_nnls = true;
solve_minimization = 1;
solve_HS = 1;
compute_Anorm = 1;

% select algorithm parameters
window_type = 'triangular';
d = 512; %(power(2, nlevel)-1)*(2*8 - 2); % amount of overlap (absolute number of pixels in the overlap)

irt_library = ['nufft_', job_id, '/'];

% path to data files prepared by Abdullah (band-interleaved)
% addpath ../data_mnras_dr
addpath /home/shared/sc004/real_data_dr
visibility_file_name = 'CYG_data_raw_ind=';

% path for data extraction
extraction_path = '/home/shared/sc004/cyg_data_dr';
% extraction_path = 'real_data_dr/';

% path to source codes
addpath ../fouRed/
addpath ../lib/
addpath ../lib/generate_data/
addpath ../lib/operators/
addpath ../lib/utils/
addpath ../lib/Proximity_operators/code/matlab/indicator/
addpath ../lib/Proximity_operators/code/matlab/multi/
addpath ../sdwt2/
addpath ../src/
addpath ../src/spmd
addpath ../src/spmd/dr/
addpath(irt_library)

% reduction parameters
% usingReduction = 1;
% usingReductionPar = 0;
% usingPrecondition = 1;
% normalize_data = 0;
% klargestpercent = 20;

%% Config parameters
Nx = 2560;
Ny = 1536;
N = Nx * Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

%% Extract / load real data (includes Fourier reduction and NNLS)
get_real_data_dr

%% Compute MAP estimator
if solve_minimization     
    Solver_fouRed_real_data
end
