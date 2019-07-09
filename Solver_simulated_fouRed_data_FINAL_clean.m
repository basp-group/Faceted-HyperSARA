%% Compute full measurement operator spectral norm
if compute_Anorm
%     if usingReduction
    if usingPrecondition
        F = afclean( @(x) HS_fouRed_forward_operator_new(x, A, At, H, Ti, Wm, aW));
        Ft = afclean( @(y) HS_fouRed_adjoint_operator_new(y, A, At, H, Ti, Wm, [Ny, Nx], aW));
    else
        F = afclean( @(x) HS_fouRed_forward_operator_new(x, A, At, H, Ti, Wm));
        Ft = afclean( @(y) HS_fouRed_adjoint_operator_new(y, A, At, H, Ti, Wm, [Ny, Nx]));
    end
%     Anorm = pow_method_op(F, Ft, [Ny Nx length(ch)]); 
    Anorm = 21174485950.1509;
%     else
%         if exist(['./simulated_data/data/Anorm.mat'], 'file')
%             load(['./simulated_data/data/Anorm.mat']);
%         else
%             F = afclean( @(x) HS_forward_operator_precond(x, Gw, A, aW));
%             Ft = afclean( @(y) HS_adjoint_operator_precond(y, Gw, At, aW, Ny, Nx));
%             Anorm = pow_method_op(F, Ft, [Ny Nx length(ch)]);
%             save(['./simulated_data/data/Anorm.mat'],'-v7.3', 'Anorm');
%         end
%     end
end

clear F Ft;
% %% Generate initial epsilons by performing imaging with NNLS on each data block separately
% if generate_eps_nnls
%     % param_nnls.im = im; % original image, used to compute the SNR
%     param_nnls.verbose = 2; % print log or not
%     param_nnls.rel_obj = 1e-5; % stopping criterion
%     param_nnls.max_iter = 1000; % max number of iterations
%     param_nnls.sol_steps = [inf]; % saves images at the given iterations
%     param_nnls.beta = 1;
%     % solve nnls per block
%     if load_data
%         load('./simulated_data/data/eps.mat');
%     else
%         for i = 1:length(ch)
%             eps_b{i} = cell(length(G{i}),1);
%             for j = 1 : length(G{i})
%                 % printf('solving for band %i\n\n',i)
%                 [~,eps_b{i}{j}] = fb_nnls_blocks(yb{i}{j}, A, At, G{i}{j}, W{i}{j}, param_nnls);
%             end
%         end
%         if save_data
%             save('./simulated_data/data/eps.mat','-v7.3', 'eps_b');
%         end
%     end
% end

%% sparsity operator definition
nlevel = 4; % wavelet level
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used
L = [2*(1:8)'; 0]; % length of the filters (0 corresponding to the 'self' basis)

if flag_algo < 2
    
    [Psi1, Psit1] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    P = length(Psi1);
    
    for k = 1 : P
        f = '@(y) HS_forward_sparsity(y,Psi1{';
        f = sprintf('%s%i},Ny,Nx);', f,k);
        Psi{k} = eval(f);
        
        ft = '@(x) HS_adjoint_sparsity(x,Psit1{';
        ft = sprintf('%s%i},%i);', ft,k,1);
        Psit{k} = eval(ft);
    end    
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
param_HSI.max_iter = 1000; % max number of iterations

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
    
%% L21 + Nuclear (facet-based version)
if solve_HS
    
    param_HSI2 = param_HSI;
    param_HSI2.nu2 = Anorm; % bound on the norm of the operator A*G
    param_HSI.reweight_alpha_ff = 0.9;
    param_HSI2.reweight_abs_of_max = 0.005;
    param_HSI2.use_reweight_steps = 1;
    param_HSI2.total_reweights = 30;
    param_HSI2.use_reweight_eps = 0;
    
    % spectral tesselation (non-overlapping)
    rg_c = domain_decomposition(Qc2, ch(end));
    cell_c_chunks = cell(Qc2, 1);
    y_spmd = cell(Qc2, 1);
    epsilon_spmd = cell(Qc2, 1);
    aW_spmd = cell(Qc2, 1);
    W_spmd = cell(Qc2, 1);
    T_spmd = cell(Qc2, 1);
    
    for i = 1:Qc2
        cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
        y_spmd{i} = yT(cell_c_chunks{i});
        epsilon_spmd{i} = epsilont(cell_c_chunks{i});
        aW_spmd{i} = aW(cell_c_chunks{i});
        W_spmd{i} = Wm(cell_c_chunks{i});
        T_spmd{i} = Ti(cell_c_chunks{i});
    end
    clear yT epsilont aW Wm Ti epsilons_t
    % old solver
%     [xsol,v0,v1,v2,g,weights0,weights1,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res] = ...
%         pdfb_LRJS_Adapt_blocks_rwNL21_par_sing_precond_new_sim
%     (yT, [Ny, Nx, length(ch)], epsilons_t, FIpsf, FIpsf_t, aW, T, W, Psi, Psit, param_HSI,X0);
    
    % new solvers
    [xsol,v0,v1,v2,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
        pdfb_LRJS_precond_NL21_sdwt2_spmd4_dr(y_spmd, [Ny, Nx], epsilon_spmd, ...
        A, At, H, aW_spmd, T_spmd, W_spmd, param_HSI2, X0, Qx, Qy, Qc2, ...
        wlt_basis, L, nlevel, cell_c_chunks, ch(end));
    
%     [xsol,v0,v1,v2,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
%     pdfb_LRJS_precond_NL21_sdwt2_spmd4_cst_overlap_weighted_dr(y_spmd, [Ny, Nx], ...
%     epsilon_spmd, A, At, H, aW_spmd, T_spmd, W_spmd, param_HSI2, X0, Qx, Qy, Qc2, ...
%     wlt_basis, L, nlevel, cell_c_chunks, ch(end), d, window_type);
    
%     [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
%     pdfb_LRJS_precond_NL21_sdwt2_spmd4_cst_overlap_weighted_dr_real(y_spmd, [Ny, Nx], ...
%     epsilon_spmd, A, At, H, aW_spmd, T_spmd, W_spmd, param_HSI2, Qx, Qy, Qc2, ...
%     wlt_basis, L, nlevel, cell_c_chunks, ch(end), d, window_type);
    
    c = size(xsol,3);
    sol = reshape(xsol(:),numel(xsol(:))/c,c);
    SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)))
    psnrh = zeros(c,1);
    for i = 1:length(ch)
        psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
    end
    SNR_average = mean(psnrh)    
    
%     mkdir('results/')
%     save(['results/results_hyperSARA_fouRed_', alg_version, '_', parallel_version, '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc), '.mat'],'-v7.3','xsol', 'sol', 'X0', 'SNR', 'SNR_average', 'res');
%     fitswrite(xsol,['results/x_hyperSARA_fouRed_', alg_version, '_', parallel_version, '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc), '.fits'])

end
