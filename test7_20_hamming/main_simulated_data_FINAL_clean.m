function main_simulated_data_FINAL(flag_algo,parallel_version,gamma,job_id)

%clear all;
%close all;
%clc;

T = 1500;
hrs = 6;

img_size = 1024;
tot = 20;
q1 = 4;
q2 = 20;

Qc = 1;
ind = 0;
Qx = q1;
Qy = q1;
Qc2 = q2;

window_type = 'hamming';
parallel_version = 'spmd4_cst_weighted';

if strcmp(parallel_version, 'spmd4_cst_weighted') || strcmp(parallel_version, 'spmd4_cst')
    nlevel = 4;
    d = (power(2, nlevel)-1)*(2*8 - 2); % assuming db8 largest wavelet filter
end


%flag_algo = 2;
%parallel_version = 'spmd4';

%gamma = 1e-4;

irt_library = ['../nufft_', job_id, '/']

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

%try
%    run('../irt/setup.m');
%end

seed = 1;
rng(seed);

generate_groundtruth = 1; % Default 60 channels cube with line emissions
generate_undersampled_cube = 0; % Default 15 channels cube with line emissions

generate_simple_image = 0; % Default the W28 image
consider_simple_image = 0; % Consider only one channel of the generated cube

generate_data = 1;

generate_uv = 1;
save_uv = 0;
load_uv = 0;

save_data = 0;
load_data = 0;

save_full_operator = 0;
load_full_operator = 0;

free_memory = 0;
free_Gw = 0;

generate_measurements = 1;
load_measurements = 0;

solve_minimization = 1;


compute_Anorm_preconditioned = 1;
save_Anorm_preconditioned = 0;

compute_Anorm = 1;
save_Anorm = 0;

compute_Anorm_ch = 0;
save_Anorm_ch = 0;

generate_eps_nnls = 0;

%% Generate or Load ground-truth
if generate_groundtruth
    file_name = ['W28_', num2str(img_size), '.fits'];
    Nx=img_size;
    Ny=img_size;
    if tot >= 60
        c = tot;
    else
        c = 60;
        generate_undersampled_cube = 1;
    end
    f = linspace(1,2,c); % frequency bandwidth from 1 to 2 GHz with 15 channels
    emission_lines = 1; % insert emission lines on top of the continuous spectra
    [x0,X0] = Generate_cube(file_name,Ny,Nx,f,emission_lines);
    ch = [1 : c]; % number of channels loaded (note that this can be one).
end

%%
if generate_undersampled_cube
    % new desired dimensions
    Nx=img_size;
    Ny=img_size;
    % undersampling factor for the channels
    unds = c/tot; % take 1/unds images
    [x0,X0,f,c] = Generate_undersampled_cube_new(x0,f,Ny,Nx,c,unds);
    ch = [1 : c]; % number of channels loaded (note that this can be one).
    Nx
    Ny
    c
end

% %% Display spectra
% figure, hold on;
% for i = 1 : 10 : size(X0,1)
%     if X0(i,end) > 0
%         plot (X0(i,:))
%     end
% end
% xlim([1 c]);

%% For monochromatic imaging
%
if generate_simple_image
    x0 = fitsread(['W28_', num2str(img_size), '.fits']);
    X0 = x0(:);
    ch = 1;
    f = 1.4;
    Nx=img_size;
    Ny=img_size;
end

%
if consider_simple_image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ch = 15; % This allows to consider one channel of the generated cube
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Generate facets
% spectral faceting (interlaved sampling)
id = interleaved_facets(Qc, c);

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

% create the complete tessellation (spectral)
if ind > 0
    x0 = x0(:,:,id{ind});
    c = size(x0,3);
    ch = [1:c];
    f = f(id{ind});
    X0 = reshape(x0,Nx*Ny,c);
end

% create the complete tessellation (spatial)
for q = 1:Q
    x_split{q} = x0(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2),:);
end

%
%if flag_algo == 0 % L11
%    param_HSI.num_workers = 34;
%elseif flag_algo == 1 % HyperSARA
%    param_HSI.num_workers = 1;
%elseif flag_algo == 2 % Faceted HyperSARA
%    param_HSI.num_workers = length(x_split)*2+1;  %%%%%%%%%%% TO BE SET BY P.-A.
%end

if flag_algo == 1 % HyperSARA - serial
    param_HSI.num_workers = 1
end

if flag_algo == 11 % HyperSARA - parfeval
    param_HSI.num_workers = 12
end

%param_HSI.num_workers


%% Generate or Load subsampling mask
if generate_data
    Generate_data_new
end

%%
sigma_noise = 0.1; %Equivalent to having inSNR = 40 dB
% Generate or Load subsampling mask
if generate_measurements
    Generate_Measurements
end

clear c;
%% Compute map estimator
if solve_minimization
    Solver_simulated_data_FINAL_clean
end

% end
%