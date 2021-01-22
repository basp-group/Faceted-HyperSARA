%% Generate initial epsilons by performing imaging with NNLS on each data block separately
if generate_eps_nnls
    % param_nnls.im = im; % original image, used to compute the SNR
    param_nnls.verbose = 2; % print log or not
    param_nnls.rel_obj = 1e-5; % stopping criterion
    param_nnls.max_iter = 1000; % max number of iterations
    param_nnls.sol_steps = [inf]; % saves images at the given iterations
    param_nnls.beta = 1;
    % solve nnls per block
    util_create_pool(32);
    parfor i = ch
        eps_b{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            % printf('solving for band %i\n\n',i)
            [~,eps_b{i}{j}] = fb_nnls_blocks(yb{i}{j}, A, At, G{i}{j}, W{i}{j}, param_nnls);
        end
    end
    save(['eps_ind=' num2str(ind) '.mat'],'-v7.3', 'eps_b');
else
    load(['eps_ind=' num2str(ind) '.mat']);
end

%%
if compute_Anorm
    %Compute full measurement operator spectral norm
    F = afclean( @(x) HS_forward_operator_precond_G(x, G, W, A, aW));
    Ft = afclean( @(y) HS_adjoint_operator_precond_G(y, G, W, At, aW, Ny, Nx));
    %disp('y')
    %size(yb{1}{1})
    %disp('G')
    %size(G{1}{1})
    %disp('W')
    %size(W{1}{1})
    %disp('aW')
    %size(aW{1}{1})
    %disp('Fx')
    %Fx = A(zeros(Ny,Nx))
    %size(FX)   
    %disp('Fx(W)')
    %size(Fx(W{1}{1}))

    Anorm = pow_method_op(F, Ft, [Ny Nx length(ch)]);
    save(['Anorm_ind=' num2str(ind) '.mat'],'-v7.3', 'Anorm');
else
    load(['Anorm_ind=' num2str(ind) '.mat']); 
end

%% sparsity operator definition
nlevel = 4; % wavelet level
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used, always put self in last position if used
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
end

% %% Splitting operator FULL
% Sp = @(x) Split_forward_operator(x,I,dims,Q);
% Spt = @(x) Split_adjoint_operator(x,I,dims,Q,Ny,Nx,length(ch));
% 
% Sp_norm = pow_method_op(Sp,Spt,[Ny Nx length(ch)]); % [P.-A.] in theory, Sp_pnorm = Pnorm (no need to compute both)

%% HSI parameter structure sent to the  HSI algorithm
param_HSI.verbose = 2; % print log or not
param_HSI.nu0 = 1; % bound on the norm of the Identity operator
param_HSI.nu1 = 1; % bound on the norm of the operator Psi
param_HSI.gamma0 = gamma0;
param_HSI.gamma = gamma;  %convergence parameter L1 (soft th parameter)
param_HSI.rel_obj = 1e-6; % stopping criterion
param_HSI.max_iter = 100000; % max number of iterations

param_HSI.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ballparam_HSI.adapt_eps_tol_out = 1.001; % tolerance outside the l2 ball
param_HSI.adapt_eps_tol_out = 1.005; % tolerance outside the l2 ball
param_HSI.adapt_eps_start = 500; % minimum num of iter before stating adjustment
param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
param_HSI.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
param_HSI.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data

param_HSI.reweight_alpha = (0.8)^alpha; 1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
param_HSI.reweight_alpha_ff = 0.8; 0.9;
param_HSI.total_reweights = 20; -1; % -1 if you don't want reweighting
param_HSI.reweight_abs_of_max = Inf; 0.005; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised

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

alpha
%% HyperSARA-sdwt2 (split L21 + nuclear norms + wavelets)

if flag_algo==2
    param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G 
    disp('Split L21 + Nuclear + wavelets')
    % spectral tesselation (non-overlapping)
    rg_c = domain_decomposition(Qc2, length(ch));
    cell_c_chunks = cell(Qc2, 1);
    y_spmd = cell(Qc2, 1);
    epsilon_spmd = cell(Qc2, 1);
    aW_spmd = cell(Qc2, 1);
    W_spmd = cell(Qc2, 1);
    G_spmd = cell(Qc2, 1);
    
    for i = 1:Qc2
        cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
        y_spmd{i} = yb(cell_c_chunks{i});
        epsilon_spmd{i} = eps_b(cell_c_chunks{i});
        aW_spmd{i} = aW(cell_c_chunks{i});
        W_spmd{i} = W(cell_c_chunks{i});
        G_spmd{i} = G(cell_c_chunks{i});
    end
    
    clear yb eps_b aW W G

    if  rw >= 0

       init_file_name = ['./results/result_HyperSARA_spmd4_cst_weighted_rd_' num2str(param_HSI.ind) '_' num2str(param_HSI.gamma0) '_' num2str(param_HSI.gamma) '_' num2str(rw) '_' num2str(alpha) '.mat']; 

    switch parallel_version

        case 'spmd4_cst_weighted'% same as spmd_sct, weight correction (apodization window in this case)
            [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
                pdfb_LRJS_precond_NL21_sdwt2_spmd4_cst_overlap_weighted_rd_abd(y_spmd, epsilon_spmd, ...
                A, At, aW_spmd, G_spmd, W_spmd, param_HSI, Qx, Qy, Qc2, ...
                wlt_basis, L, nlevel, cell_c_chunks, length(ch), d, window_type,init_file_name);

        otherwise
            error('Unknown parallelisation option.')
    end

    else

    switch parallel_version

        case 'spmd4_cst_weighted'% same as spmd_sct, weight correction (apodization window in this case)
            [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
                pdfb_LRJS_precond_NL21_sdwt2_spmd4_cst_overlap_weighted_rd_abd(y_spmd, epsilon_spmd, ...
                A, At, aW_spmd, G_spmd, W_spmd, param_HSI, Qx, Qy, Qc2, ...
                wlt_basis, L, nlevel, cell_c_chunks, length(ch), d, window_type, '');
            
        otherwise
            error('Unknown parallelisation option.')
    end

    end
    
    Time_iter_average = mean(end_iter)
    
    mkdir('results/')
    save(['results/results_hyperSARA_', parallel_version, '_ind=', num2str(ind), '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc2), '_gamma=', num2str(gamma),'.mat'],'-v7.3','xsol', 'param', 'epsilon', 't', 'rel_fval', 'nuclear', 'l21', 'norm_res', 'res','end_iter');

end