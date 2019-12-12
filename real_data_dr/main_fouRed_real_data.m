clear all
close all
clc

%% Paramters
chInd = 32; %1:32
reduction_version = 2; % 1 for F\Phi^t, 2 for G^t
realdatablocks = 2;
enable_klargestpercent = false;
fouRed_gamma = 3;

solve_minimization = 1;
gamma = 1e-6;
algo_version = 2; % 1 for facet HyperSARA with DR, 2 for SARA with DR

%% Extract / load real data (includes Fourier reduction and NNLS)
if reduction_version
    for j = 1:length(chInd)
        func_nnls_dr_real_data(chInd(j), reduction_version, realdatablocks, enable_klargestpercent, fouRed_gamma)
    end
    fprintf('Reduction finished\n');
end

%% Compute MAP estimator
if solve_minimization 
    if reduction_version
        func_solver_fouRed_real_data(gamma, chInd, reduction_version, algo_version, realdatablocks, fouRed_gamma)
    else
        func_solver_real_data(gamma, chInd)
    end
end
