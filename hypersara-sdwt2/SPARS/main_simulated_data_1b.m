function main_simulated_data(flag_algo, tot, chunk_width, overlap)

% clear all;
% close all;
% clc;

Qx = 3; % additional parameters for the faceting
Qy = 3;
%         
% num_workers = 12;
% tot = 14;
% chunk_width = 5;
% step = 3;
%param_HSI.num_workers = num_workers;

step = chunk_width - overlap;
chunks = [];
for  i = 1 :step: tot-chunk_width+1
    chunks = [chunks; i i+chunk_width-1];
end
diff = tot - chunks(end,2);
if diff <= 0.4*chunk_width
    chunks(end,2) = tot;
else
    chunks = [chunks; chunks(end,2) tot];
end

chunks

if flag_algo == 0
param_HSI.num_workers = 34;
elseif flag_algo == 1
param_HSI.num_workers = 1;
elseif flag_algo == 3
param_HSI.num_workers = size(chunks,1)*2+1;
else
param_HSI.num_workers = size(chunks,1)+1;
end
param_HSI.num_workers

% clear all;
% close all;
% clc;

addpath ../simulated_data/Audrey/
addpath ../simulated_data/data/
addpath ../simulated_data/
addpath ../lib/
addpath ../lib/CubeHelix/
addpath ../lib/irt/
addpath ../lib/sdwt2
addpath ../lib/Proximity_operators/code/matlab/indicator/
addpath ../lib/Proximity_operators/code/matlab/multi/
addpath ../alg

try
    run('../lib/irt/setup.m');
end

generate_groundtruth = 1; % Default 60 channels cube with line emissions
load_groundtruth = 0;

generate_undersampled_cube = 0; % Default 15 channels cube with line emissions

generate_simple_image = 0; % Default the W28 image
consider_simple_image = 0; % Consider only one channel of the generated cube

generate_data = 1;

gaussian_UV = 0;
realistic_UV = 1;

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

% choose the solver you want to use
solve_HS = 0; % wide-band solver (rwLRJS) applied on all the channels
solve_minimization_UQ_HS = 0;

load_xmap = 1;
save_xmap = 0;

load_mask = 1;
save_mask = 0;

compute_Anorm_preconditioned = 0;
load_Anorm_preconditioned = 0;
save_Anorm_preconditioned = 0;

compute_Anorm = 0;
load_Anorm = 0;
save_Anorm = 0;

solve_1B = 1; % single-channel solver (SARA) applied for each channel separately
solve_minimization_UQ_1B = 0;

compute_Anorm_ch = 1;
save_Anorm_ch = 0;
load_Anorm_ch = 0;

generate_eps_nnls = 0;
load_temp_result = 0;

%% Generate or Load ground-truth
if generate_groundtruth
    file_name = '../simulated_data/data/W28_256.fits';
    Nx=256;
    Ny=256;
    c = tot;
    f = linspace(1,2,c); % frequency bandwidth from 1 to 2 GHz with 15 channels
    emission_lines = 0; % insert emission lines on top of the continuous spectra
    [x0,X0] = Generate_cube(file_name,Ny,Nx,f,emission_lines);
    ch = [1 : c]; % number of channels loaded (note that this can be one).
elseif load_groundtruth
    load('../simulated_data/data/data_60_HI_final3.mat','M','N','f','c','xf','XF');
    Nx=M;
    Ny=N;
    x0 = xf;
    X0 = XF;
    clear M N XF xf;
    ch = [1 : c]; % number of channels loaded (note that this can be one).
end

%%
if generate_undersampled_cube
    % new desired dimensions
    Nx=4*1024;
    Ny=4*1024;
    % undersampling factor for the channels
    unds = 1; % take 1/unds images
    [x0,X0,f,c] = Generate_undersampled_cube(x0,f,Ny,Nx,c,unds);
    ch = [1 : c]; % number of channels loaded (note that this can be one).
end

%%
if generate_simple_image
    x0 = fitsread('W28_256.fits');
    X0 = x0(:);
    ch = 1;
    f = 1.4;
    Nx=256;
    Ny=256;
end

%%
if consider_simple_image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ch = 15; % This allows to consider one channel of the generated cube
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%
percentage = 1; [0.3 0.05 0.1 0.02 0.01 0.2]; % sampling rate
for k = 1 : length(percentage)
    
    seed = 1;
    rng(seed);
    
    %diary(['diary_HI_60_' num2str(percentage(k))]);
    
    %% Generate or Load subsampling mask
    if generate_data
        Generate_data
    end
    
    %%
    sigma_noise = 0.1;
    input_snr = 100; % noise level (on the measurements)
    num_tests = 1;
    % Generate or Load subsampling mask
    if generate_measurements
        Generate_Measurements
        if save_data
            save(['../simulated_data/data/y_60_HI_finalNEW_' num2str(percentage(k)) '.mat'],'-v7.3', 'y0_t', 'y_t', 'y0b_t', 'yb_t', 'epsilonT', 'epsilonT1', 'sigma_noise');
        end
    elseif load_measurements
        load(['../simulated_data/data/y_60_HI_finalNEW_' num2str(percentage(k)) '.mat'],'y0b_t','yb_t');
        disp('data loaded successfully ;)')
        
    end
    clear c;
    
    %% Compute map estimator
    if solve_minimization
        Solver_simulated_data_splitting_FINAL_sdwt2
    end
    
    %% Compute uncertainty quantification
    
    if solve_minimization_UQ_HS
        
        Solver_UQ_simulated_data_HS_crop
        
    end
end

%