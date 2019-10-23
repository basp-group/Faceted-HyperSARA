clear all
close all
clc

%% Paramters
chInd = 32;
reduction_version = 2; % 1 for F\Phi^t, 2 for G^t

solve_minimization = 1;
gamma = 1e-6;
algo_version = 2; % 1 for HyperSARA, 2 for SARA

%% Extract / load real data (includes Fourier reduction and NNLS)
if reduction_version
    nnls_dr_reduced_data(chInd, reduction_version)
%     dr_real_data(chInd, reduction_version)
    fprintf('Reduction finished\n');
end

%% Compute MAP estimator
if solve_minimization 
    if reduction_version
        script_half_Solver_fouRed_real_data(gamma, chInd, reduction_version, algo_version)
    else
        script_Solver_real_data(gamma, chInd)
    end
end
