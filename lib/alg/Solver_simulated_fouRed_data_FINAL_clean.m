fprintf('Qx=%d, Qy=%d, Qc2=%d\n', Qx, Qy, Qc2);
d = 256;
rw = -1;
window_type = 'triangular';
realdatablocks = length(H{1});

cell_c_chunks = cell(Qc2, 1);
nchannel_per_worker = zeros(Qc2, 1);
rg_c = domain_decomposition(Qc2, length(ch)); % to be modified % only for data reduction
for i = 1:Qc2
    cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
%     cell_c_chunks{i} = [rg_c1(i, 1):rg_c1(i, 2), length(ch)-rg_c2(i, 2)+1:length(ch)-rg_c2(i, 1)+1];   % only for data reduction
    nchannel_per_worker(i) = numel(cell_c_chunks{i});
end

%% Load data
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

Ap = Composite();
Atp = Composite();
Hp = Composite();
Wp = Composite();
yTp = Composite();
Tp = Composite();
aWp = Composite();
Wmp = Composite();
epsilonp = Composite();
l2_upper_boundp = Composite();

for i = 1:Qc2
    Ap{i+Q} = A;
    Atp{i+Q} = At;
    Hp{i+Q} = H(cell_c_chunks{i});
    Wp{i+Q} = W(cell_c_chunks{i});
    yTp{i+Q} = yT(cell_c_chunks{i});
    Tp{i+Q} = T(cell_c_chunks{i});
    aWp{i+Q} = aW(cell_c_chunks{i});
    Wmp{i+Q} = Wm(cell_c_chunks{i});
    epsilonp{i+Q} = epsilons(cell_c_chunks{i});
    l2_upper_boundp{i+Q} = epsilons(cell_c_chunks{i});
end
clear A At H Wm yT T aW Wm epsilons

%% Compute full measurement operator spectral norm
Anormfile = './data/Anorm.mat';

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
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8'}; % wavelet basis to be used
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

%% L21 + Nuclear (facet-based version)
if solve_HS
    disp('Split L21 + Nuclear + wavelets')    
    %% HSI parameter structure sent to the  HSI algorithm
    param_HSI.verbose = 2; % print log or not
    param_HSI.nu0 = 1; % bound on the norm of the Identity operator
    param_HSI.nu1 = 1; % bound on the norm of the operator Psi
    param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
    param_HSI.gamma0 = 1;
    param_HSI.gamma = 1e-4;  %convergence parameter L1 (soft th parameter)
    param_HSI.rel_obj = 1e-10; % stopping criterion
    param_HSI.max_iter = 1000; % max number of iterations

    param_HSI.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
    param_HSI.adapt_eps_start = 300; % minimum num of iter before stating adjustment
    param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
    param_HSI.adapt_eps_tol_out = 1.01; % tolerance outside the l2 ball
    param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
    param_HSI.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
    param_HSI.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data

    param_HSI.reweight_alpha = (0.8)^20; %1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
    param_HSI.reweight_alpha_ff = 0.63;
    param_HSI.total_reweights = 5; % -1 if you don't want reweighting
    param_HSI.reweight_abs_of_max = Inf; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised

    param_HSI.use_reweight_steps = 0; % reweighting by fixed steps
    param_HSI.reweight_step_size = 300; % reweighting step size
    param_HSI.reweight_steps = [5000: param_HSI.reweight_step_size :10000];
    param_HSI.step_flag = 1;

    param_HSI.use_reweight_eps = 1; % reweighting w.r.t the relative change of the solution
    param_HSI.reweight_max_reweight_itr = param_HSI.max_iter - param_HSI.reweight_step_size;
    param_HSI.reweight_rel_obj = 5e-4; % criterion for performing reweighting
    param_HSI.reweight_min_steps_rel_obj = 300; % min num of iter between reweights

    param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid
    param_HSI.precondition = usingPrecondition;
    param_HSI.rw_tol = 300;
    param_HSI.ind = 1:1;
    
%     param_HSI.initsol = xsol;
    
    reweight_step_count = 0;
    initfilename = ['./results/facetHyperSARA_dr_co_w_real_' ...
                num2str(param_HSI.ind(1)), '_', num2str(param_HSI.ind(end)), '_' num2str(param_HSI.gamma) '_' num2str(reweight_step_count) '.mat'];
        
    % spectral tesselation (non-overlapping)
%     epsilon_spmd = cell(Qc2, 1);
       
%     for i = 1:Qc2
%         if adapt_eps_flag
%             epsilon_spmd{i} = epsilon1(cell_c_chunks{i});
%         else
%         epsilon_spmd{i} = epsilon(cell_c_chunks{i});
%         end
%         l2_upper_bound{i} = epsilon(cell_c_chunks{i});
%     end
    
%     param_HSI.l2_upper_bound = l2_upper_bound;
    
%     clear epsilon epsilon1
    
    if  rw >= 0 
        load(['results/result_HyperSARA_spmd4_cst_weighted_rd_' num2str(param_HSI.gamma) '_' num2str(rw) '.mat']);
        
%         if adapt_eps_flag
%             epsilon_spmd = epsilon1;
%         else
%             epsilon_spmd = epsilon;
%         end
        param_HSI.init_xsol = param.init_xsol;
        param_HSI.init_g = param.init_g;
        param_HSI.init_v0 = param.init_v0;
        param_HSI.init_v1 = param.init_v1;
        param_HSI.init_weights0 = param.init_weights0;
        param_HSI.init_weights1 = param.init_weights1;
        param_HSI.init_v2 = param.init_v2;
        param_HSI.init_proj = param.init_proj;
        param_HSI.init_t_block = param.init_t_block;
        param_HSI.init_t_start = param.init_t_start+1;
        param_HSI.reweight_alpha = param.reweight_alpha;
        param_HSI.init_reweight_step_count = param.init_reweight_step_count;
        param_HSI.init_reweight_last_iter_step = param.init_reweight_last_iter_step;
        param_HSI.reweight_steps = (param_HSI.init_t_start+param_HSI.reweight_step_size: param_HSI.reweight_step_size :param_HSI.max_iter+(2*param_HSI.reweight_step_size));
        param_HSI.step_flag = 0;
    end
    
    % solvers
    mkdir('results/')
    
    [xsol,param_HSI,t,rel_fval,nuclear,l21,end_iter,snr,snr_avg] = ...
        facetHyperSARA_DR_SIMU_precond_v3(yTp, epsilonp, Ap, Atp, Hp, Wp, aWp, ...
        Tp, Wmp, param_HSI, X0, Qx, Qy, Qc2, wlt_basis, L, nlevel, cell_c_chunks, ...
        length(ch), d, window_type, initfilename, reduction_version, realdatablocks,...
        fouRed_gamma, Ny, Nx, l2_upper_boundp);
    
    save(['results/results_facethyperSARA_fouRed_ch', num2str(ch(1)), '_', num2str(ch(end)), '_', num2str(flag_algo), '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), ...
        '_Qc=', num2str(Qc2), '_gamma=', num2str(gamma), '_gamma0=', num2str(gamma0), '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_adpteps', num2str(adapt_eps_flag),'.mat'], '-v7.3', ...
        'xsol', 'param_HSI', 't', 'rel_fval', 'nuclear', 'l21', 'end_iter','snr','snr_avg');
end


% %% sparsity operator definition
% nlevel = 4; % wavelet level
% wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used
% L = [2*(1:8)'; 0]; % length of the filters (0 corresponding to the 'self' basis)
% 
% if flag_algo < 2
%     
%     [Psi1, Psit1] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
%     P = length(Psi1);
%     
%     for k = 1 : P
%         f = '@(y) HS_forward_sparsity(y,Psi1{';
%         f = sprintf('%s%i},Ny,Nx);', f,k);
%         Psi{k} = eval(f);
%         
%         ft = '@(x) HS_adjoint_sparsity(x,Psit1{';
%         ft = sprintf('%s%i},%i);', ft,k,1);
%         Psit{k} = eval(ft);
%     end    
% end
% 
% %%
% % [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
% % Psi = @(y) HS_forward_sparsity(y,Psiw,Ny,Nx);
% % Psit = @(x) HS_adjoint_sparsity(x,Psitw,length(wlt_basis));
% 
% % Pnorm = pow_method_op(Psit,Psi,[Ny Nx c]);
% 
% %%
% % Anorm1 = pow_method_op(@(x) Gw{1}*A(x), @(x) At(Gw{1}' * x), [Ny Nx c]);
% 
% %% HSI parameter structure sent to the  HSI algorithm
% param_HSI.verbose = 2; % print log or not
% param_HSI.nu0 = 1; % bound on the norm of the Identity operator
% param_HSI.nu1 = 1; % bound on the norm of the operator Psi
% param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
% param_HSI.gamma0 = 1;
% param_HSI.gamma = 1e-2;  %convergence parameter L1 (soft th parameter)
% param_HSI.rel_obj = 1e-5; % stopping criterion
% param_HSI.max_iter = 1000; % max number of iterations
% 
% param_HSI.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
% param_HSI.adapt_eps_start = 500; % minimum num of iter before stating adjustment
% param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
% param_HSI.adapt_eps_tol_out = 1.001; % tolerance outside the l2 ball
% param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
% param_HSI.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
% param_HSI.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data
% 
% param_HSI.reweight_alpha = 1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
% param_HSI.reweight_alpha_ff = 0.5;
% param_HSI.total_reweights = 5; % -1 if you don't want reweighting
% param_HSI.reweight_abs_of_max = 1; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised
% 
% param_HSI.use_reweight_steps = 1; % reweighting by fixed steps
% param_HSI.reweight_step_size = 300; % reweighting step size
% param_HSI.reweight_steps = [5000: param_HSI.reweight_step_size :10000];
% param_HSI.step_flag = 1;
% 
% param_HSI.use_reweight_eps = 0; % reweighting w.r.t the relative change of the solution
% param_HSI.reweight_max_reweight_itr = param_HSI.max_iter - param_HSI.reweight_step_size;
% param_HSI.reweight_rel_obj = 5e-4; % criterion for performing reweighting
% param_HSI.reweight_min_steps_rel_obj = 300; % min num of iter between reweights
% 
% param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
% param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
% param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid
% 
% param_HSI.precondition = usingPrecondition;
% 
% %%
%     
% %% L21 + Nuclear (facet-based version)
% if solve_HS
%     
%     param_HSI2 = param_HSI;
%     param_HSI2.nu2 = Anorm; % bound on the norm of the operator A*G
%     param_HSI.reweight_alpha_ff = 0.9;
%     param_HSI2.reweight_abs_of_max = 0.005;
%     param_HSI2.use_reweight_steps = 1;
%     param_HSI2.total_reweights = 30;
%     param_HSI2.use_reweight_eps = 0;
%     
%     % spectral tesselation (non-overlapping)
%     rg_c = domain_decomposition(Qc2, ch(numel(ch)));
%     cell_c_chunks = cell(Qc2, 1);
%     y_spmd = cell(Qc2, 1);
%     epsilon_spmd = cell(Qc2, 1);
%     aW_spmd = cell(Qc2, 1);
%     W_spmd = cell(Qc2, 1);
%     T_spmd = cell(Qc2, 1);
%     
%     for i = 1:Qc2
%         cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
%         y_spmd{i} = yT(cell_c_chunks{i});
%         epsilon_spmd{i} = epsilont(cell_c_chunks{i});
%         aW_spmd{i} = aW(cell_c_chunks{i});
%         W_spmd{i} = Wm(cell_c_chunks{i});
%         T_spmd{i} = Ti(cell_c_chunks{i});
%     end
%     clear yT epsilont aW Wm Ti epsilons_t
%     % old solver
% %     [xsol,v0,v1,v2,g,weights0,weights1,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res] = ...
% %         pdfb_LRJS_Adapt_blocks_rwNL21_par_sing_precond_new_sim
% %     (yT, [Ny, Nx, length(ch)], epsilons_t, FIpsf, FIpsf_t, aW, T, W, Psi, Psit, param_HSI,X0);
%     
%     % new solvers
%     [xsol,v0,v1,v2,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
%         pdfb_LRJS_precond_NL21_sdwt2_spmd4_dr(y_spmd, [Ny, Nx], epsilon_spmd, ...
%         A, At, H, aW_spmd, T_spmd, W_spmd, param_HSI2, X0, Qx, Qy, Qc2, ...
%         wlt_basis, L, nlevel, cell_c_chunks, ch(end));
%     
%     % Serial Fourier-reduced hyperSARA
%     
% %     [xsol,v0,v1,v2,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
% %     pdfb_LRJS_precond_NL21_sdwt2_spmd4_cst_overlap_weighted_dr(y_spmd, [Ny, Nx], ...
% %     epsilon_spmd, A, At, H, aW_spmd, T_spmd, W_spmd, param_HSI2, X0, Qx, Qy, Qc2, ...
% %     wlt_basis, L, nlevel, cell_c_chunks, ch(end), d, window_type);
%     
% %     [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
% %     pdfb_LRJS_precond_NL21_sdwt2_spmd4_cst_overlap_weighted_dr_real(y_spmd, [Ny, Nx], ...
% %     epsilon_spmd, A, At, H, aW_spmd, T_spmd, W_spmd, param_HSI2, Qx, Qy, Qc2, ...
% %     wlt_basis, L, nlevel, cell_c_chunks, ch(end), d, window_type);
%     
%     c = size(xsol,3);
%     sol = reshape(xsol(:),numel(xsol(:))/c,c);
%     SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)))
%     psnrh = zeros(c,1);
%     for i = 1:length(ch)
%         psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
%     end
%     SNR_average = mean(psnrh)    
%     
% %     mkdir('results/')
% %     save(['results/results_hyperSARA_fouRed_', alg_version, '_', parallel_version, '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc2), '.mat'],'-v7.3','xsol', 'sol', 'X0', 'SNR', 'SNR_average', 'res');
% %     fitswrite(xsol,['results/x_hyperSARA_fouRed_', alg_version, '_', parallel_version, '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc2), '.fits'])
% 
% end
