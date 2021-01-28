function func_solver_fouRed_simulated_data(yT, cubeSize, A, At, H, W, T, ...
    aW, Wm, epsilon, gamma0, gamma, ch, subInd, reduction_version, ...
    algo_version, realdatablocks, fouRed_gamma, fouRed_type,...
    adapt_eps_flag, Qx, Qy, Qc2, wterm, levelG, levelC, init_file_name, rw_alpha, rw_tot, lowRes)

Ny = cubeSize(1);
Nx = cubeSize(2);
nChannels = cubeSize(3);

cell_c_chunks = cell(Qc2, 1);
nchannel_per_worker = zeros(Qc2, 1);
rg_c = domain_decomposition(Qc2, nChannels); % to be modified % only for data reduction
for i = 1:Qc2
    cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
%     cell_c_chunks{i} = [rg_c1(i, 1):rg_c1(i, 2), nChannels-rg_c2(i, 2)+1:nChannels-rg_c2(i, 1)+1];   % only for data reduction
    nchannel_per_worker(i) = numel(cell_c_chunks{i});
end

delete(gcp('nocreate'));
Q = Qx*Qy; 
numworkers = Q + Qc2;
cirrus_cluster = parcluster('local');
cirrus_cluster.NumWorkers = numworkers;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
if cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads > ncores
    exit(1);
end
parpool(cirrus_cluster, numworkers, 'IdleTimeout', Inf);

spmd
    if labindex > Q
        Ap = A;
        Atp = At;
        chunk = cell_c_chunks{labindex - Q};
        ch_len = length(chunk);
        Hp = cell(ch_len,1);
        Wp = cell(ch_len,1);
        yTp = cell(ch_len,1);
        Tp = cell(ch_len,1);
        aWp = cell(ch_len,1);
        Wmp = cell(ch_len,1);
        epsilonp = cell(ch_len,1);
        l2_upper_boundp = cell(ch_len,1);
        for i = 1: ch_len
            ch_ind = chunk(i);
            fprintf('\nChannel number: %d\n', ch(ch_ind))
            Hp{i,1} = H{ch_ind,1};
            Wp{i,1} = W{ch_ind,1};
            yTp{i,1} = yT{ch_ind,1};
            Tp{i,1} = T{ch_ind,1};
            aWp{i,1} = aW{ch_ind,1};
            Wmp{i,1} = Wm{ch_ind,1};
            epsilonp{i,1} = epsilon{ch_ind,1};
            l2_upper_boundp{i,1} = epsilon{ch_ind,1};                   
            fprintf('\nDR file of channel number %d has been read\n', ch(ch_ind))
        end
    end
end

clear H W yT T aW Wm epsilon

%% Compute full measurement operator spectral norm
if wterm
    Anormfile = ['./data/Anorm_dr_ch', num2str(ch(1)), '_', num2str(ch(end)), '_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '.mat'];
else
    Anormfile = ['./data/Anorm_dr_ch', num2str(ch(1)), '_', num2str(ch(end)), '_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.mat'];
end
if isfile(Anormfile)
    compute_Anorm = false;
else
    compute_Anorm = true;
end

if compute_Anorm
    fprintf('\nCompute the operator norm: \n')
    Anorm = pow_method_op_composite(Hp, Wp, Ap, Atp, Tp, aWp, Q, Qc2, nchannel_per_worker, [Ny Nx]);
    save(Anormfile,'-v7.3', 'Anorm');
else
    fprintf('\nLoad the operator norm file: %s\n', Anormfile)
    load(Anormfile);
end
% delete(gcp('nocreate'));
fprintf('\nThe operator norm: %e\n', Anorm)

clear F Ft;

%% Sparsity operator definition
nlevel = 4; % wavelet level
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used
L = [2*(1:8)'; 0]; % length of the filters (0 corresponding to the 'self' basis)

if flag_algo < 2
    
    [Psi1, Psit1] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    P = length(Psi1);

    for k = 1 : P
        f = '@(x_wave) HS_forward_sparsity(x_wave,Psi1{';
        f = sprintf('%s%i},Ny,Nx);', f,k);
        Psi{k} = eval(f);

        b(k) = size(Psit1{k}(zeros(Ny,Nx,1)),1);

        ft = ['@(x) HS_adjoint_sparsity(x,Psit1{' num2str(k) '},b(' num2str(k) '));'];
        Psit{k} = eval(ft);
    end

    %% Full sparsity operator
    [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    bb = size(Psitw(zeros(Ny,Nx,1)),1);

    Psi_full = @(x_wave) HS_forward_sparsity(x_wave,Psiw,Ny,Nx);
    Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,bb);
    
elseif flag_algo == 2
    [Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);    
end
%%
% [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
% Psi = @(y) HS_forward_sparsity(y,Psiw,Ny,Nx);
% Psit = @(x) HS_adjoint_sparsity(x,Psitw,length(wlt_basis));

% Pnorm = pow_method_op(Psit,Psi,[Ny Nx c]);

%%
% Anorm1 = pow_method_op(@(x) Gw{1}*A(x), @(x) At(Gw{1}' * x), [Ny Nx c]);

%% HSI parameter structure sent to the  HSI algorithm
param_HSI.verbose = 2; % print log or not
param_HSI.nu0 = 1; % bound on the norm of the Identity operator
param_HSI.nu1 = 1; % bound on the norm of the operator Psi
param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
param_HSI.gamma0 = 1;
param_HSI.gamma = 1e-2;  %convergence parameter L1 (soft th parameter)
param_HSI.rel_obj = 1e-5; % stopping criterion
param_HSI.max_iter = 500; % max number of iterations

param_HSI.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
param_HSI.adapt_eps_start = 500; % minimum num of iter before stating adjustment
param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
param_HSI.adapt_eps_tol_out = 1.001; % tolerance outside the l2 ball
param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
param_HSI.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
param_HSI.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data

param_HSI.reweight_alpha = 1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
param_HSI.reweight_alpha_ff = 0.5;
param_HSI.total_reweights = 5; % -1 if you don't want reweighting
param_HSI.reweight_abs_of_max = 1; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised

param_HSI.use_reweight_steps = 1; % reweighting by fixed steps
param_HSI.reweight_step_size = 300; % reweighting step size
param_HSI.reweight_steps = [5000: param_HSI.reweight_step_size :10000];
param_HSI.step_flag = 1;

param_HSI.use_reweight_eps = 0; % reweighting w.r.t the relative change of the solution
param_HSI.reweight_max_reweight_itr = param_HSI.max_iter - param_HSI.reweight_step_size;
param_HSI.reweight_rel_obj = 5e-4; % criterion for performing reweighting
param_HSI.reweight_min_steps_rel_obj = 300; % min num of iter between reweights

param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid

param_HSI.precondition = usingPrecondition;

%%
for q = 1 : length(sigma_noise) * num_tests %number of tests x number of InSNRs
    
%% L21 + Nuclear
%     [xsol,v0,v1,v2,g,weights0,weights1,t_block,reweight_alpha,epsilon,iterh,rel_fval,nuclear,l21,norm_res,res] = pdfb_LRJS_Adapt_blocks_rwNL21_par_sing_new_sim(yT, [Ny, Nx, length(ch)], epsilons_t, FIpsf, FIpsf_t, T, W, Psi, Psit, param_HSI, X0);
[xsol,v0,v1,v2,g,weights0,weights1,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res] = pdfb_LRJS_Adapt_blocks_rwNL21_par_sing_precond_new_sim(yT, [Ny, Nx, length(ch)], epsilons_t, FIpsf, FIpsf_t, aW, T, W, Psi, Psit, param_HSI,X0);
[xsol,param,t,rel_fval,nuclear,l21,end_iter,snr,snr_avg] = ...
    facetHyperSARA_fouRed_simu(yp, epsilonp, Ap, Atp, Hp, Wp, pUp, Tp, Wmp, ...
    param, X0, Qx, Qy, K, wavelet, L, nlevel, c_chunks, c, d, window_type, init_file_name, ...
    reduction_version, realdatablocks, fouRed_gamma, M, N, l2_upper_bound);
c = size(xsol,3);
sol = reshape(xsol(:),numel(xsol(:))/c,c);
SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)))
psnrh = zeros(c,1);
for i = 1:length(ch)
    psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
end
SNR_average = mean(psnrh)

end