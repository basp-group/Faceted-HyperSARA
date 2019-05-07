clear all
close all
clc


flag_algo = 2; tot = 20; Qx = 2; Qy = 2; Qc = 2; ind = 1; img_size = 512; %2048;
Qc2 = 2;
T = 200;
hrs = 6;
parallel_version = 'spmd4';

addpath fouRed/
addpath hypersara-clean/lib/
addpath hypersara-clean/lib/generate_data/
addpath hypersara-clean/lib/operators/
addpath hypersara-clean/lib/CubeHelix/
addpath nufft
addpath sdwt2/
addpath src/
addpath data/
addpath hypersara-clean/lib/Proximity_operators/code/matlab/indicator/
addpath hypersara-clean/lib/Proximity_operators/code/matlab/multi/

% addpath ../hypersara/simulated_data/data/
% addpath ../hypersara/simulated_data/
% addpath ../hypersara/alg/

% try
%     run('../irt/setup.m');
% end

usingReduction = 1;
usingReductionPar = 0;
usingPrecondition = 1;
normalize_data = 0;
klargestpercent = 20;


generate_groundtruth = 1; % Default 60 channels cube with line emissions
load_groundtruth = 0;
generate_undersampled_cube = 1;

simple_test_1B = 0;
generate_simple_image = 0; % Default the W28 image
consider_simple_image = 0; % Consider only one channel of the generated cube

generate_data = 1;

gaussian_UV = 1;
realistic_UV = 0;

generate_uv = 1;
save_uv = 0;
load_uv = 0;

save_data = 0;
load_data = 0;

save_full_operator = 0;
load_full_operator = 0;

free_memory = 1;
free_Gw = 0;

generate_measurements = 1;
load_measurements = 0;

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
if generate_groundtruth
    file_name = ['W28_', num2str(img_size), '.fits'];
    Nx=img_size;
    Ny=img_size;
    if tot > 60
        c = tot;
    else
        c = 60;
        generate_undersampled_cube = 1;
    end
    f = linspace(1,2,c); % frequency bandwidth from 1 to 2 GHz with 15 channels
    emission_lines = 1; % insert emission lines on top of the continuous spectra
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

%%
if consider_simple_image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ch = 15; % This allows to consider one channel of the generated cube
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Generate facets
% spectral faceting (interlaved sampling)
id = interleaved_facets(Qc, c);

% create the complete tessellation (spectral)
if ind > 0
    x0 = x0(:,:,id{ind});
    c = size(x0,3);
    ch = [1:c];
    f = f(id{ind});
    X0 = reshape(x0,Nx*Ny,c);
end

%
if flag_algo == 0 % L11
    param_HSI.num_workers = 34;
elseif flag_algo == 1 % HyperSARA
    param_HSI.num_workers = 1;
elseif flag_algo == 2 % Faceted HyperSARA
    param_HSI.num_workers = Qx*Qy*2+1;  %%%%%%%%%%% TO BE SET BY P.-A.
end
param_HSI.num_workers


%% Generate or Load subsampling mask
if generate_data
    Generate_data_new
end


percentage = 0.3; %[0.3 0.05 0.1 0.02 0.01 0.2]; % sampling rate

%% TO BE MODIFIED FROM HERE

%% Generate or Load subsampling mask
if usingReduction
    Generate_fouRed_data
else
    Generate_data
end

%%
sigma_noise = 0.1;
%     input_snr = 40; % noise level (on the measurements)
% Generate or Load subsampling mask
if usingReduction
    Generate_fouRed_Measurements
else 
    if generate_measurements
        Generate_Measurements
%             if save_data
%                 save(['../simulated_data/data/y_60_HI_finalNEW_' num2str(percentage(k)) '.mat'],'-v7.3', 'y0_t', 'y_t', 'y0b_t', 'yb_t', 'epsilonT', 'epsilonT1', 'sigma_noise');
%             end
%         elseif load_measurements
%             load(['../simulated_data/data/y_60_HI_finalNEW_' num2str(percentage(k)) '.mat'],'y0b_t','yb_t');
%             disp('data loaded successfully ;)')

    end
end
clear c;
    
%% Fourier reduction
if usingReduction
    Fourier_reduction_ch
end

%% Compute MAP estimator
if solve_minimization     
    if usingReduction
        Solver_simulated_fouRed_data_FINAL_clean
    end
end
