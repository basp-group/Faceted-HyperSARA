%% Compute full measurement operator spectral norm
if compute_Anorm
    F = afclean( @(x) HS_forward_operator_precond(x, Gw, A, aW));
    Ft = afclean( @(y) HS_adjoint_operator_precond(y, Gw, At, aW, Ny, Nx));
    Anorm = pow_method_op(F, Ft, [Ny Nx length(ch)]);
    if save_data
        save('./simulated_data/data/Anorm.mat','-v7.3', 'Anorm');
    end
elseif load_data
    load('./simulated_data/data/Anorm.mat');
end

%% Compute measurement operator spectral norm for each channel individually
if compute_Anorm_ch
    for i = ch
        Anorm_ch(i) = pow_method_op(@(x) sqrt(cell2mat(aW{i})) .* (Gw{i}*A(x)), @(x) At(Gw{i}' * (sqrt(cell2mat(aW{i})) .* x)), [Ny Nx 1]);
    end
    if save_data
        save('./simulated_data/data/Anorm_ch.mat','-v7.3', 'Anorm_ch');
    end
elseif load_data
    load('./simulated_data/data/Anorm_ch.mat');
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
        save('./simulated_data/data/eps.mat','-v7.3', 'eps_b');
    end
elseif load_data
    load('./simulated_data/data/eps.mat');
end

%% sparsity operator definition
nlevel = 4; % wavelet level
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used
% wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8'}; % wavelet basis to be used

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

%% Splitting operator
tot = length(ch);
chunk_width = 1/2 * tot;
step = 1/4 * tot;

chunks = [];
for  i = 1 :step: tot-chunk_width+1
    chunks = [chunks; i i+chunk_width-1];
end

Sp = @(x) Split_forward_operator(x,chunks);
Spt = @(x) Split_adjoint_operator2(x,chunks,Ny,Nx,length(ch));

Sp_norm = pow_method_op(Sp,Spt,[Ny Nx length(ch)]);

%% Average sparsity operator 
[Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
Psi_full = @(y) HS_forward_sparsity(y,Psiw,Ny,Nx);
Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,length(wlt_basis));

Sp_psit = @(x) Split_sparsity_forward_operator(x,Sp,Psit_full);
Spt_psi = @(x) Split_sparsity_adjoint_operator(x,Spt,Psi_full);

Pnorm = pow_method_op(Sp_psit, Spt_psi,[Ny Nx length(ch)]);

%%
% Anorm1 = pow_method_op(@(x) Gw{1}*A(x), @(x) At(Gw{1}' * x), [Ny Nx c]);

%% HSI parameter structure sent to the  HSI algorithm
param_HSI.verbose = 2; % print log or not
param_HSI.nu0 = Sp_norm; % bound on the norm of the Identity operator
param_HSI.nu1 = Pnorm; % bound on the norm of the operator Psi
param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
param_HSI.gamma0 = 1;
param_HSI.gamma = 1e-2;  %convergence parameter L1 (soft th parameter)
param_HSI.rel_obj = 5e-4; % stopping criterion
param_HSI.max_iter = 2000; % max number of iterations

param_HSI.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
param_HSI.adapt_eps_start = 200; % minimum num of iter before stating adjustment
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

parameters_option = 1;

if parameters_option == 1
    
    normalisation_type = 'No_normalised_operator';
    
    param_HSI.ll = 1;
    param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
    
    y_t1 = y_t;
    epsilons_t1 = epsilons_t;
    
elseif parameters_option == 2
    
    normalisation_type = 'normalised_operator';
    %
    param_HSI.ll = 1/sqrt(Anorm);
    for i = ch
        for j = 1 : length(y_t{1}{i})
            y_t1{1}{i}{j} = param_HSI.ll * y_t{1}{i}{j};
            epsilons_t1{1}{i}{j} = param_HSI.ll * epsilons_t{1}{i}{j};
        end
    end
    param_HSI.nu2 = Anorm * (param_HSI.ll)^2;
    
end

%%
for q = 1 : length(input_snr) * num_tests %number of tests x number of InSNRs
    
%% L21 + Nuclear
if solve_HS
    [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,iterh,rel_fval,nuclear,l21,norm_res,res] = pdfb_LRJS_Adapt_blocks_rwNL21_par_precond_new_sim_NL21_split(y_t{q}, epsilons_t{q}, A, At, aW, G, W, Sp, Spt, Psi, Psit, param_HSI, X0);
    if save_data
        save('./simulated_data/data/result.mat','-v7.3','xsol', 'v0', 'v1', 'v2', 'g', 'weights0', 'weights1', 'proj', 't_block','reweight_alpha', 'epsilon', 'iterh', 'rel_fval', 'nuclear', 'l21', 'norm_res', 'res');
    end
    
    c = size(xsol,3);
    sol = reshape(xsol(:),numel(xsol(:))/c,c);
    SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)))
    psnrh = zeros(c,1);
    for i = ch
        psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
    end
    SNR_average = mean(psnrh)
    
end

%% rwL11 Adaptive Blocks
if solve_1B
    for i = ch
        i
        param_HSI.nu2 = Anorm_ch(i); % bound on the norm of the operator A*G
        [xsol,v1,v2,g,weights1,proj,t_block,reweight_alpha,epsilon,iterh,rel_fval,l11,norm_res,res] = pdfb_L11_Adapt_blocks_rw_par_precond_new({y_t{q}{i}},{epsilons_t{q}{i}}, A, At, {aW{i}}, {G{i}}, {W{i}}, Psi, Psit, param_HSI, i);
        if save_data
            save(['./simulated_data/data/result_b' num2str(i) '.mat'],'-v7.3', 'xsol', 'v1', 'v2', 'g', 'weights1', 'proj', 't_block','reweight_alpha', 'epsilon', 'iterh', 'rel_fval', 'l11', 'norm_res', 'res');
        end
    end
end

end

%%
save('result.mat','-v7.3','xsol', 'sol', 'X0', 'SNR', 'SNR_average', 'res');

fitswrite(xsol,'xsol.fits')
fitswrite(x0,'x0.fits')

c