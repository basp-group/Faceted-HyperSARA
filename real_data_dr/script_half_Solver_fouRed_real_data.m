function script_half_Solver_fouRed_real_data(gamma, ch, reduction_version, algo_version)

addpath ../fouRed
addpath ../lib/
addpath ../lib/operators/
addpath ../lib/nufft/
addpath ../lib/utils/
addpath ../lib/CubeHelix/
addpath ../lib/Proximity_operators/code/matlab/indicator/
addpath ../lib/Proximity_operators/code/matlab/multi/
addpath ../sdwt2/
addpath ../src/
addpath ../src/spmd/
addpath ../src/spmd/dr/
addpath ../src/spmd/weighted/

fprintf("gamma=%e\n", gamma)
disp(ch)

compute_Anorm = true;
usingPrecondition = false;
rw = -1;
window_type = 'triangular';

Qx = 5;
Qy = 3;
Qc2 = 1;    % to see later

d = 512;
flag_algo = 2;
parallel_version = 'spmd4_cst';
bool_weights = true; % for the spmd4_new version (50% overlap version)

param_real_data.image_size_Nx = 2560; % 2560;
param_real_data.image_size_Ny = 1536; % 1536;
nChannels = length(ch); % total number of "virtual" channels (i.e., after
% concatenation) for the real dataset considered
nBlocks = 2;        % number of data blocks (needs to be known beforehand,
% quite restrictive here), change l.70 accordingly
% klargestpercent = 20;
FT2 = @(x) fftshift(fft2(ifftshift(x)));

%% Config parameters
Nx = param_real_data.image_size_Nx;
Ny = param_real_data.image_size_Ny;
N = Nx * Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

[A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);

%% Load data
for i = 1:nChannels 
    ch(i)
    tmp = load(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_H=', num2str(ch(i)), '.mat']);
%     tmp = load(['/home/basphw/mjiang/Data/mjiang/real_data_dr/CYG_old_yT=', num2str(i), '.mat'], 'yT');
%     tmp = load(['./CYG_yT=', num2str(ch(i)), '.mat'], 'yT');
    H{i,1} = tmp.H{1,1};
    tmp = load(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_DR_th3=', num2str(ch(i)), '.mat']);
%     tmp = load(['/home/basphw/mjiang/Data/mjiang/real_data_dr/CYG_old_DR=', num2str(i), '.mat']);
%     tmp = load(['./CYG_DR=', num2str(ch(i)), '.mat']);
    yT{i,1} = tmp.yT{1,1};
    T{i,1} = tmp.T{1,1};
    aW{i,1} = tmp.aW{1,1};
    Wm{i,1} = tmp.Wm{1,1};
%     tmp = load(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_gam01_epsilon=', num2str(ch(i)), '.mat'], 'epsilon');
%     tmp = load(['/home/basphw/mjiang/Data/mjiang/real_data_dr/CYG_old_epsilon=', num2str(i), '.mat'], 'epsilon');
%     tmp = load(['./CYG_epsilon=', num2str(ch(i)), '.mat'], 'epsilon');
    epsilon{i,1} = tmp.epsilon{1,1};
end
% yT{1}(1) = [];
% H{1}(1) = [];
% T{1}(1) = [];
% aW{1}(1) = [];
% Wm{1}(1) = [];
% epsilon{1}(1) = [];

clear tmp

%% Compute full measurement operator spectral norm
if compute_Anorm
    if reduction_version == 1
        F = afclean( @(x) HS_fouRed_forward_operator_new(x, A, At, H, T, Wm, aW));
        Ft = afclean( @(y) HS_fouRed_adjoint_operator_new(y, A, At, H, T, Wm, [Ny, Nx], aW));
    elseif reduction_version == 2
        F = afclean( @(x) HS_operatorGtPhi(x, H{1}, A, T{1}, Wm{1}));
        Ft = afclean( @(y) HS_operatorGtPhi_t(y, H{1}, At, T{1}, Wm{1}, [Ny, Nx], [oy * Ny, ox * Nx]));
    end
    Anorm = pow_method_op(F, Ft, [Ny Nx nChannels]);    
    save(['Anorm_dr_0prec_per_', num2str(ch(1)), '_', num2str(ch(end)),'_th3.mat'],'-v7.3', 'Anorm');
else
    load(['Anorm_dr_0prec_per_', num2str(ch(1)), '_', num2str(ch(end)),'.mat']);
end

Anorm

clear F Ft;

%% Sparsity operator definition
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
    
    % Full sparsity operator
    [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    bb = size(Psitw(zeros(Ny,Nx,1)),1);
    
    Psi_full = @(x_wave) HS_forward_sparsity(x_wave,Psiw,Ny,Nx);
    Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,bb);
    % Pnorm = pow_method_op(Psit,Psi,[Ny Nx c]);
    
elseif flag_algo == 2
    [Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);    
end

%% L21 + Nuclear (facet-based version)
if algo_version == 1
    disp('Split L21 + Nuclear + wavelets')    
    %% HSI parameter structure sent to the  HSI algorithm
    param_HSI.verbose = 2; % print log or not
    param_HSI.nu0 = 1; % bound on the norm of the Identity operator
    param_HSI.nu1 = 1; % bound on the norm of the operator Psi
    param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
    param_HSI.gamma0 = 1;
    param_HSI.gamma = gamma;  %convergence parameter L1 (soft th parameter)
    param_HSI.rel_obj = 1e-10; % stopping criterion
    param_HSI.max_iter = 100000; % max number of iterations

    param_HSI.use_adapt_eps = 1; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
    param_HSI.adapt_eps_start = 300; % minimum num of iter before stating adjustment
    param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
    param_HSI.adapt_eps_tol_out = 1.001; % tolerance outside the l2 ball
    param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
    param_HSI.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
    param_HSI.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data

    param_HSI.reweight_alpha = 1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
    param_HSI.reweight_alpha_ff = 0.8;
    param_HSI.total_reweights = 40; % -1 if you don't want reweighting
    param_HSI.reweight_abs_of_max = 1; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised

    param_HSI.use_reweight_steps = 1; % reweighting by fixed steps
    param_HSI.reweight_step_size = 300; % reweighting step size
    param_HSI.reweight_steps = [5000: param_HSI.reweight_step_size :10000];
    param_HSI.step_flag = 1;

    param_HSI.use_reweight_eps = 1; % reweighting w.r.t the relative change of the solution
    param_HSI.reweight_max_reweight_itr = param_HSI.max_iter - param_HSI.reweight_step_size;
    param_HSI.reweight_rel_obj = 1e-4; % criterion for performing reweighting
    param_HSI.reweight_min_steps_rel_obj = 300; % min num of iter between reweights

    param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid
    param_HSI.precondition = usingPrecondition;

    
    % spectral tesselation (non-overlapping)
    rg_c = domain_decomposition(Qc2, nChannels); % to be modified
    cell_c_chunks = cell(Qc2, 1);
    y_spmd = cell(Qc2, 1);
    epsilon_spmd = cell(Qc2, 1);
    aW_spmd = cell(Qc2, 1);
    W_spmd = cell(Qc2, 1);
    T_spmd = cell(Qc2, 1);
    
    for i = 1:Qc2
        cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
        y_spmd{i} = yT(cell_c_chunks{i});
        epsilon_spmd{i} = epsilon(cell_c_chunks{i});
        aW_spmd{i} = aW(cell_c_chunks{i});
        W_spmd{i} = Wm(cell_c_chunks{i});
        T_spmd{i} = T(cell_c_chunks{i});
        H_spmd{i} = H(cell_c_chunks{i});
    end
    clear yT epsilon aW Wm T epsilon
    
    if  rw >= 0 
        load(['results/result_HyperSARA_spmd4_cst_weighted_rd_' num2str(param_HSI.gamma) '_' num2str(rw) '.mat']);
        
        epsilon_spmd = epsilon;
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
    [xsol,param_HSI,epsilon,t,rel_fval,nuclear,l21,norm_res_out,res,end_iter] = ...
        facetHyperSARA_cst_overlap_weighted_dr_real_data(y_spmd, [Ny, Nx], ...
        epsilon_spmd, A, At, H_spmd, aW_spmd, T_spmd, W_spmd, param_HSI, Qx, Qy, Qc2, ...
        wlt_basis, L, nlevel, cell_c_chunks, nChannels, d, window_type);
    
    save(['results/results_hyperSARA_fouRed_ch', num2str(ch(1)), '_', num2str(ch(end)),'_', algo_version, '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), ...
        '_Qc=', num2str(Qc2), '_gamma=', num2str(gamma),'.mat'], '-v7.3', ...
        'xsol', 'param_HSI', 'epsilon', 't', 'rel_fval', 'nuclear', 'l21', ...
        'norm_res', 'res', 'end_iter');

elseif algo_version == 2
    disp('Split L1 + wavelets')
    % %% PDFB parameter structure sent to the algorithm
    param_pdfb.verbose = 1; % print log or not
    param_pdfb.nu1 = 1; % bound on the norm of the operator Psi
    param_pdfb.nu2 = Anorm; % bound on the norm of the operator A*G
    param_pdfb.gamma = gamma; % convergence parameter L1 (soft th parameter)
    param_pdfb.tau = 0.49; % forward descent step size
    param_pdfb.rel_obj = 1e-6; % stopping criterion
    param_pdfb.max_iter = 10000; % max number of iterations
    param_pdfb.lambda0 = 1; % relaxation step for primal update
    param_pdfb.lambda1 = 1; % relaxation step for L1 dual update
    param_pdfb.lambda2 = 1; % relaxation step for L2 dual update
    param_pdfb.omega1 = 1;
    param_pdfb.omega2 = 1;
    param_pdfb.sol_steps = [inf]; % saves images at the given iterations

    param_pdfb.use_proj_elipse_fb = logical(usingPrecondition);
    param_pdfb.elipse_proj_max_iter = 200;
    param_pdfb.elipse_proj_min_iter = 1;
    param_pdfb.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid
    
    param_pdfb.use_adapt_eps = 1; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
    param_pdfb.adapt_eps_start = 300; % minimum num of iter before stating adjustment
    param_pdfb.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
    param_pdfb.adapt_eps_tol_out = 1.001; % tolerance outside the l2 ball
    param_pdfb.adapt_eps_steps = 100; % min num of iter between consecutive updates
    param_pdfb.adapt_eps_rel_obj = 1e-3; % bound on the relative change of the solution
    param_pdfb.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data

    param_pdfb.reweight_alpha = 1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
    param_pdfb.reweight_alpha_ff = 0.8;
    param_pdfb.total_reweights = 40; % -1 if you don't want reweighting
    param_pdfb.reweight_abs_of_max = 1; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised

    param_pdfb.use_reweight_steps = 1; % reweighting by fixed steps
    param_pdfb.reweight_step_size = 300; % reweighting step size
    param_pdfb.reweight_steps = [2000: param_pdfb.reweight_step_size :10000];
    param_pdfb.step_flag = 1;

    param_pdfb.use_reweight_eps = 1; % reweighting w.r.t the relative change of the solution
    param_pdfb.reweight_max_reweight_itr = param_pdfb.max_iter - param_pdfb.reweight_step_size;
    param_pdfb.reweight_rel_obj = 1e-4; % criterion for performing reweighting
    param_pdfb.reweight_min_steps_rel_obj = 300; % min num of iter between reweights
    
    % solvers
    mkdir('results/')
    [xsol, t_block, epsilon, t, rel_fval, L1_v, L2_v, L2_vp, norm2, res, end_iter] = ...
        pdfb_bpcon_DR_adapt_eps_precond(yT{1}, [Ny, Nx], epsilon{1}, A, At, H{1}, aW{1}, T{1}, Wm{1}, ... 
        Psi, Psit, Psiw, Psitw, param_pdfb, reduction_version);
    
    save(['results/results_fouRed_th3_ch', num2str(ch(1)), '_', num2str(ch(end)),'_', algo_version, '_gamma=', num2str(gamma),'.mat'], '-v7.3', ...
        'xsol', 'epsilon', 't', 'rel_fval', 'norm2', 'res', 'end_iter');
    
end
end
