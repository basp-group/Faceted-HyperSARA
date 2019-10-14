clear all
close all
clc

%% Paramters
chInd = 1;
reduction_version = 1; % 1 for F\Phi^t, 2 for G^t

solve_minimization = 1;
gamma = 1e-6;
algo_version = 2; % 1 for HyperSARA, 2 for SARA

%% Extract / load real data (includes Fourier reduction and NNLS)
dr_real_data(chInd, reduction_version)
fprintf('Reduction finished\n');

%% Compute MAP estimator
if solve_minimization     
    script_half_Solver_fouRed_real_data(gamma, [1:1], reduction_version, algo_version)
end
