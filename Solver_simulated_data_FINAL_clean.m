if flag_algo == 0 % solve L11 problem
    
    % Compute measurement operator spectral norm for each channel individually
    for i = ch
        Anorm_ch(i) = pow_method_op(@(x) sqrt(cell2mat(aW{i})) .* (Gw{i}*A(x)), @(x) At(Gw{i}' * (sqrt(cell2mat(aW{i})) .* x)), [Ny Nx 1]);
    end
    
else
    if compute_Anorm_preconditioned
        % Compute full measurement operator spectral norm
        F = afclean( @(x) HS_forward_operator_precond(x, Gw, A, aW));
        Ft = afclean( @(y) HS_adjoint_operator_precond(y, Gw, At, aW, Ny, Nx));
        Anorm = pow_method_op(F, Ft, [Ny Nx length(ch)]);
        if save_data
            mkdir('simulated_data/data/')
            save('simulated_data/data/Anorm.mat','-v7.3', 'Anorm');
        end
    else
        load('simulated_data/data/Anorm.mat','-v7.3', 'Anorm');
    end
end

if exist('Gw','var')
    clear Gw;
end

%% Generate initial epsilons by performing imaging with NNLS on each data block separately
if generate_eps_nnls
    % param_nnls.im = im; % original image, used to compute the SNR
    param_nnls.verbose = 2; % print log or not
    param_nnls.rel_obj = 1e-5; % stopping criterion
    param_nnls.max_iter = 1000; % max number of iterations
    param_nnls.sol_steps = [inf]; % saves images at the given iterations
    param_nnls.beta = 1;
    % solve nnls per block
    for i = ch
        eps_b{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            % printf('solving for band %i\n\n',i)
            [~,eps_b{i}{j}] = fb_nnls_blocks(yb{i}{j}, A, At, G{i}{j}, W{i}{j}, param_nnls);
        end
    end
    if save_data
        mkdir('simulated_data/data')
        save('simulated_data/data/eps.mat','-v7.3', 'eps_b');
    end
elseif load_data
    load('simulated_data/data/eps.mat');
end

%% sparsity operator definition
nlevel = 4; % wavelet level
%wlt_basis = {'db8'};
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used
% wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8'}; % wavelet basis to be used
L = [2*(1:8)'; 0]; % length of the filters

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

%% Splitting operator FULL
Sp = @(x) Split_forward_operator(x,I,dims,Q);
Spt = @(x) Split_adjoint_operator(x,I,dims,Q,Ny,Nx,length(ch));

Sp_norm = pow_method_op(Sp,Spt,[Ny Nx length(ch)]); % [P.-A.] in theory, Sp_pnorm = Pnorm (no need to compute both)

%% HSI parameter structure sent to the  HSI algorithm
param_HSI.verbose = 2; % print log or not
param_HSI.nu0 = 1; % bound on the norm of the Identity operator
param_HSI.nu1 = 1; % bound on the norm of the operator Psi
param_HSI.gamma0 = 1;
param_HSI.gamma = 1e-2;  %convergence parameter L1 (soft th parameter)
param_HSI.rel_obj = 1e-5; % stopping criterion
param_HSI.max_iter = 10000; % max number of iterations

param_HSI.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
param_HSI.adapt_eps_start = 200; % minimum num of iter before stating adjustment
param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
param_HSI.adapt_eps_tol_out = 1.001; % tolerance outside the l2 ball
param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
param_HSI.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
param_HSI.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data

param_HSI.reweight_alpha = 1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
param_HSI.reweight_alpha_ff = 0.9;
param_HSI.total_reweights = 30; -1; % -1 if you don't want reweighting
param_HSI.reweight_abs_of_max = 0.005; 1; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised

param_HSI.use_reweight_steps = 1; % reweighting by fixed steps
param_HSI.reweight_step_size = 300; % reweighting step size
param_HSI.reweight_steps = [5000: param_HSI.reweight_step_size :10000];
param_HSI.step_flag = 1;

param_HSI.use_reweight_eps = 0; % reweighting w.r.t the relative change of the solution
param_HSI.reweight_max_reweight_itr = param_HSI.max_iter - param_HSI.reweight_step_size;
param_HSI.reweight_rel_obj = 1e-4; 5e-4; % criterion for performing reweighting
param_HSI.reweight_min_steps_rel_obj = 300; % min num of iter between reweights

param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid

%% HyperSARA
if flag_algo == 1
    param_HSI2 = param_HSI;
    param_HSI2.nu2 = Anorm; % bound on the norm of the operator A*G
    param_HSI.reweight_alpha_ff = 0.9;
    param_HSI2.reweight_abs_of_max = 0.005;
    param_HSI2.use_reweight_steps = 1;
    param_HSI2.total_reweights = 30;
    param_HSI2.use_reweight_eps = 0;
    
    disp('-----------------------------------------')
    disp('HyperSARA')
    disp('-----------------------------------------')
    
%     [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,iterh,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
%         pdfb_LRJS_Adapt_blocks_rwNL21_par_precond_new_sim(y, epsilons, A, At, aW, G, W, Psi, Psit, param_HSI2, X0);
    
    % test
    rg_c = domain_decomposition(Qc, ch(end));
    cell_c_chunks = cell(Qc, 1);
    y_spmd = cell(Qc, 1);
    epsilon_spmd = cell(Qc, 1);
    aW_spmd = cell(Qc, 1);
    W_spmd = cell(Qc, 1);
    G_spmd = cell(Qc, 1);
    
    for i = 1:Qc
        cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
        y_spmd{i} = y(cell_c_chunks{i});
        epsilon_spmd{i} = epsilons(cell_c_chunks{i});
        aW_spmd{i} = aW(cell_c_chunks{i});
        W_spmd{i} = W(cell_c_chunks{i});
        G_spmd{i} = G(cell_c_chunks{i});
    end
    [xsol,v0,v1,v2,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
        pdfb_LRJS_precond_NL21_sdwt2_spmd_serial_SARA(y_spmd, epsilon_spmd, A, At, aW_spmd, G_spmd, W_spmd, param_HSI2, X0, Qc, wlt_basis, nlevel, cell_c_chunks, ch(end));
    % - end test
    
    c = size(xsol,3);
    sol = reshape(xsol(:),numel(xsol(:))/c,c);
    SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)))
    psnrh = zeros(c,1);
    for i = ch
        psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
    end
    SNR_average = mean(psnrh)
    
    Time_iter_average = mean(end_iter)
    
    mkdir('results')
    save(['results/result_HyperSARA_' num2str(T) '_' num2str(hrs) '.mat'],'-v7.3','xsol', 'sol', 'X0', 'SNR', 'SNR_average', 'res','end_iter');
    fitswrite(xsol,'results/xsol_HyperSARA.fits')
    fitswrite(x0,'results/x0.fits')
end

%% SARA
if flag_algo == 0
    param_HSI2 = param_HSI;
    param_HSI2.use_reweight_steps = 1;
    param_HSI2.total_reweights = 30;
    util_create_pool(param_HSI.num_workers);
    maxNumCompThreads(param_HSI.num_workers);
    parfor i = 1 : length(ch)
        i
        %param_HSI.nu2 = Anorm_ch(i); % bound on the norm of the operator A*G
        [xsol(:,:,i),v1,v2,g,weights1,proj,t_block,reweight_alpha,epsilon,iterh,rel_fval,l11,norm_res,res(:,:,i),end_iter{i}] = ...
            pdfb_L11_Adapt_blocks_rw_precond_new_sim({y{i}},{epsilons{i}}, A, At, {aW{i}}, {G{i}}, {W{i}}, Psi_full, Psit_full, param_HSI2, X0(:,i), Anorm_ch(i));
        
    end
    
    c = size(xsol,3);
    sol = reshape(xsol(:),numel(xsol(:))/c,c);
    SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)))
    psnrh = zeros(c,1);
    for i = ch
        psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
        Time_iter = mean(end_iter{i});
    end
    SNR_average = mean(psnrh)
    
    Time_iter_average = mean(Time_iter)
    
    mkdir('results')
    save(['results/result_L11_' num2str(T) '_' num2str(hrs) '.mat'],'-v7.3','xsol', 'sol', 'X0', 'SNR', 'SNR_average', 'res','end_iter', 'Time_iter_average');
    fitswrite(xsol,'results/xsol_L11.fits')
    fitswrite(x0,'results/x0.fits')
end

%% HyperSARA-sdwt2 (split L21 + nuclear norms + wavelets)

if flag_algo==2
    
    param_HSI2 = param_HSI;
    param_HSI2.nu2 = Anorm; % bound on the norm of the operator A*G
    param_HSI.reweight_alpha_ff = 0.9;
    param_HSI2.reweight_abs_of_max = 0.005;
    param_HSI2.use_reweight_steps = 1;
    param_HSI2.total_reweights = 30;
    param_HSI2.use_reweight_eps = 0;
    
    %copy-paste my current codes, to be cleansed before going to cirrus
     %% 
    % [P.-A.] only portion to be revised---------')
    disp('Split L21 + Nuclear + wavelets')
    disp('-----------------------------------------')
%     [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,iterh,rel_fval,nuclear,l21,norm_res,res] = ...
%     pdfb_LRJS_Adapt_blocks_rwNL21_par_precond_sim_NL21_split_sdwt2(y_t{q}, epsilons_t{q}, A, At, aW, G, W, Sp, Spt, param_HSI, X0, Qx, Qy, wlt_basis, L, nlevel);   

    % [P.-A.]
%     [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
%     pdfb_LRJS_precond_NL21_sdwt2(y_t{q}, epsilons_t{q}, A, At, aW, G, W, Sp, Spt, param_HSI, X0, Qx, Qy, num_chunk, wlt_basis, L, nlevel, c_chunks, chunks, Psit_full);

    switch parallel_version
%         case 'parfeval'
%             [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
%             pdfb_LRJS_precond_NL21_sdwt2_parfeval(y, epsilons, A, At, aW, G, W, Sp, Spt, param_HSI, X0, Qx, Qy, num_chunk, wlt_basis, L, nlevel, c_chunks, Psit_full);
        case 'parfeval2'
            [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
            pdfb_LRJS_precond_NL21_sdwt2_parfeval2(y, epsilons, A, At, aW, G, W, param_HSI2, X0, Qx, Qy, wlt_basis, L, nlevel, Psit_full);
        case 'spmd3'
            [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
            pdfb_LRJS_precond_NL21_sdwt2_spmd3(y, epsilons, A, At, aW, G, W, param_HSI2, X0, Qx, Qy, wlt_basis, L, nlevel);
        case 'spmd4'
            % spectral tesselation (non-overlapping)
            rg_c = domain_decomposition(Qc, ch(end));
            cell_c_chunks = cell(Qc, 1);
            y_spmd = cell(Qc, 1);
            epsilon_spmd = cell(Qc, 1);
            aW_spmd = cell(Qc, 1);
            W_spmd = cell(Qc, 1);
            G_spmd = cell(Qc, 1);
            
            for i = 1:Qc
                cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
                y_spmd{i} = y(cell_c_chunks{i});
                epsilon_spmd{i} = epsilons(cell_c_chunks{i});
                aW_spmd{i} = aW(cell_c_chunks{i});
                W_spmd{i} = W(cell_c_chunks{i});
                G_spmd{i} = G(cell_c_chunks{i});
            end

            [xsol,v0,v1,v2,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
                pdfb_LRJS_precond_NL21_sdwt2_spmd4(y_spmd, epsilon_spmd, ...
                A, At, aW_spmd, G_spmd, W_spmd, param_HSI2, X0, Qx, Qy, Qc, ...
                wlt_basis, L, nlevel, cell_c_chunks, ch(end));
            
%             [xsol,v0,v1,v2,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
%                 pdfb_LRJS_precond_NL21_sdwt2_spmd5(y_spmd, epsilon_spmd, ...
%                 A, At, aW_spmd, G_spmd, W_spmd, param_HSI2, X0, Qx, Qy, Qc, ...
%                 wlt_basis, L, nlevel, cell_c_chunks, ch(end)); % gather image and data on the same nodes (extra communications compared to spmd4 for reweigthing and monitoring variables)
%             
%             [xsol,v0,v1,v2,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
%                 pdfb_LRJS_precond_NL21_sdwt2_spmd6(y_spmd, epsilon_spmd, ...
%                 A, At, aW_spmd, G_spmd, W_spmd, param_HSI2, X0, Qx, Qy, Qc, ...
%                 wlt_basis, L, nlevel, cell_c_chunks, ch(end));
        otherwise
            error('Unknown parallelisation option.')
    end
    
    c = size(xsol,3);
    sol = reshape(xsol(:),numel(xsol(:))/c,c);
    SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)))
    psnrh = zeros(c,1);
    for i = ch
        psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
    end
    SNR_average = mean(psnrh)
    Time_iter_average = mean(end_iter)
    
    mkdir('results/')
    save(['results/results_hyperSARA_', alg_version, '_', parallel_version, '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc), '.mat'],'-v7.3','xsol', 'sol', 'X0', 'SNR', 'SNR_average', 'res');
    fitswrite(xsol,['results/x_hyperSARA_', alg_version, '_', parallel_version, '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc), '.fits'])

end
