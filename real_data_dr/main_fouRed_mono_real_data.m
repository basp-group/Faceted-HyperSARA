clear all
close all
clc

%% Paramters
chInd = 32;
reduction_version = 2; % 1 for F\Phi^t, 2 for G^t
realdatablocks = 2;
enable_klargestpercent = false;
fouRed_gamma = 3;

solve_minimization = 1;
gamma = 1e-6;
algo_version = 2; % 1 for HyperSARA, 2 for SARA
realdatablocks = 2;

%% Extract / load real data (includes Fourier reduction and NNLS)
if reduction_version
    nnls_dr_reduced_data(chInd, reduction_version, realdatablocks, enable_klargestpercent, fouRed_gamma)
    fprintf('Reduction finished\n');
end

%% Compute MAP estimator
if solve_minimization 
    if reduction_version
        script_half_Solver_fouRed_real_data(gamma, chInd, reduction_version, algo_version, realdatablocks, fouRed_gamma)
    else
        script_Solver_real_data(gamma, chInd)
    end
end
