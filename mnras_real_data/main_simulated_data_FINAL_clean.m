function main_simulated_data_FINAL_clean(ind,gamma0,gamma,rw,alpha,job_id)

% clear all;
% close all;
% clc;

% ind = 5;
% gamma = 1e-5;

param_HSI.alpha = alpha;

rw_tol = 5000;
param_HSI.rw_tol = rw_tol;

%ind = 6;
param_HSI.ind = ind;

Qx = 5;
Qy = 3;
Qc2 = 15;

flag_algo = 2;

window_type = 'triangular';
parallel_version = 'spmd4_cst_weighted';

if strcmp(parallel_version, 'spmd4_cst_weighted') || strcmp(parallel_version, 'spmd4_cst')
    nlevel = 4;
    d = 512; %(power(2, nlevel)-1)*(2*8 - 2); % assuming db8 largest wavelet filter
end

irt_library = ['../nufft_', job_id, '/']
%irt_library = ['../nufft/']


addpath ../hypersara-clean/lib/
addpath ../hypersara-clean/lib/generate_data/
addpath ../hypersara-clean/lib/operators/
addpath ../hypersara-clean/lib/CubeHelix/
addpath(irt_library)
addpath ../sdwt2/
addpath ../src/
addpath ../data/
addpath ../hypersara-clean/lib/Proximity_operators/code/matlab/indicator/
addpath ../hypersara-clean/lib/Proximity_operators/code/matlab/multi/
%addpath ../hypersara-clean/alg


extract_raw_data = 0;
compute_Anorm = 0;
generate_real_data = 1;
generate_eps_nnls = 0;

save_data = 0;
save_full_operator = 0;
free_memory = 1;

solve_minimization = 1;
load_temp_result = 0;

%% Config parameters
if compute_Anorm || ~generate_real_data
    Nx = 2560;
    Ny = 1536;
    N = Nx * Ny;
    ox = 2; % oversampling factors for nufft
    oy = 2; % oversampling factors for nufft
    Kx = 8; % number of neighbours for nufft
    Ky = 8; % number of neighbours for nufft
end

%% Generate or Load real data
if generate_real_data
    extract_real_data;
else
    load('CYG_data.mat');
    ch = [1 : size(yb,2)]
    [A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);
end

%% Generate facets
% spectral faceting (interlaved sampling)
%id = interleaved_facets(Qc, c);

% spatial faceting
Q = Qx*Qy;
rg_y = domain_decomposition(Qy, Ny);
rg_x = domain_decomposition(Qx, Nx);

% create starting index of the spatial facets (I) and associated dimensions
% (dims). Beware: I starts at (0, 0)
I = zeros(Q, 2);
dims = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        I(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1];
        dims(q, :) = [rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
    end
end

%% Compute map estimator
if solve_minimization
    Solver_simulated_data_FINAL_clean
end

end
