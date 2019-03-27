clear all
close all
clc

% addpath simulated_data/data/
% addpath simulated_data/
% addpath lib/
% addpath lib/irt/
% addpath alg
addpath ../fouRed/
addpath ../hypersara/lib/
addpath ../hypersara/lib/irt/
addpath ../hypersara/simulated_data/data/
addpath ../hypersara/simulated_data/
addpath ../hypersara/alg/

try
    run('./lib/irt/setup.m');
end

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

param_HSI.num_workers = 24;

chunk_width = 2;
step = 1;

chunks = [];
for  i = 1 :step: 15-chunk_width+1
    chunks = [chunks; i i+chunk_width-1];
end

chunks

%% Generate or Load ground-truth
if generate_groundtruth
    file_name = '../hypersara/simulated_data/data/W28_256.fits';
    Nx=256;
    Ny=256;
    c = 61;
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
    Nx=256;
    Ny=256;
    % undersampling factor for the channels
    unds = 4; % take 1/unds images
    [x0,X0,f,c] = Generate_undersampled_cube(x0,f,Ny,Nx,c-1,unds);
    ch = [1 : c]; % number of channels loaded (note that this can be one).
end

%%
if consider_simple_image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ch = 15; % This allows to consider one channel of the generated cube
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%
percentage = 0.3; %[0.3 0.05 0.1 0.02 0.01 0.2]; % sampling rate

for k = 1 : length(percentage)
    
    seed = 1;
    rng(seed);
    
    %diary(['diary_HI_60_' num2str(percentage(k))]);
    
    %% Generate or Load subsampling mask
    if usingReduction
        Generate_fouRed_data
    else
        Generate_data
    end
    
    %%
    sigma_noise = 0.1;
%     input_snr = 40; % noise level (on the measurements)
    num_tests = 1;
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
    
    %% Compute map estimator
    if solve_minimization     
        if usingReduction
            Solver_fouRed_simulated_data_FINAL
        else
            Solver_simulated_data_splitting_FINAL
%             Solver_simulated_data
        end
    end
    
    %% Compute uncertainty quantification
    
    if solve_minimization_UQ_HS
        
        Solver_UQ_simulated_data_HS_crop
        
    end
end

