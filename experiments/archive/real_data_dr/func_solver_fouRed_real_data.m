function func_solver_fouRed_real_data(datadir, gamma0, gamma, ch, subInd, reduction_version, algo_version, realdatablocks, fouRed_gamma, fouRed_type, adapt_eps_flag, wterm, levelG, levelC, init_file_name, rw_alpha, rw_tot, lowRes)
% Global imaging solver for DR real data
%-------------------------------------------------------------------------%
% Input:
% > datadir: directory of data, string
% > name: filename for savings, string
% > Qx: number of facets along dimension x [1]
% > Qy: number of facets along dimension y [1]
% > Qc2: number of datab computing processes [1]
% > gamma0: parameter for nuclear norm, [1]
% > gamma1: parameter for l2,1 norm, [1]
% > ch: vector of channel indices [L]
% > subInd: vector of interlaced subproblem indices [S]
% > reduction_version: reduction version, [1]
%       1: old one (not used any more)
%       2: current one
% > algo_version: version of algorithm, [1]
%       1. hyperspectral version
%       2. monochromatic version
% > realdatablocks: number of data blocks, [1]
% > fouRed_gamma: level of reduction, [1]
% > fouRed_type: type of reduction, supplementary information for
% fouRed_gamma, [1]
%       1: remove smallest singular values based on "fouRed_gamma" percentage
%       2: remove smallest singular values based on "fouRed_gamma"-sigma
% > adapt_eps_flag: flag of adaptive epsilon strategy, [1]
% > wterm: whether consider w-term, [1]
% > levelG: thresholding level to reduce size when consider w-term, [1]
% > levelC: thresholding level to reduce size when consider w-term, [1]
% > init_file_name: file name to initialize/resume running, string
% > rw_alpha: coefficient for reweighting formula, [1]
% > rw_tot: total reweighting times, [1]
% > lowRes: flag of using low resolution image, used for tests, [1]

if fouRed_type == 1
    typeStr = 'perc';
elseif fouRed_type == 2
    typeStr = 'th';
end

addpath ../fouRed;
addpath ../lib/;
addpath ../lib/operators/;
addpath ../lib/measurement-operator/nufft/;
addpath ../lib/utils/;
addpath ../lib/faceted-wavelet-transform/src;
addpath ../src_mnras/;
addpath ../src_mnras/spmd/;
addpath ../src_mnras/spmd/dr/;
addpath ../src_mnras/spmd/weighted/;

fprintf("Nuclear norm parameter=%e\n", gamma0);
fprintf("Regularization parameter=%e\n", gamma);
fprintf('Channel number %d\n', ch);
fprintf('Index number %d\n', subInd);
fprintf('Reduction version %d\n', reduction_version);
fprintf('Algorithm version: %d\n', algo_version);
fprintf('Data blocks: %d\n', realdatablocks);
if fouRed_type == 1
    fprintf('Reduction level: remove %f percentile\n', fouRed_gamma);
elseif fouRed_type == 2
    fprintf('Reduction level: keep %f sigma\n', fouRed_gamma);
end
fprintf('Adapt epsilon: %d\n', adapt_eps_flag);
fprintf('Reweight alpha: %f, total reweights: %d\n', rw_alpha, rw_tot);
fprintf('Low resolution: %d\n', lowRes);

usingPrecondition = true;
rw = -1;
window_type = 'triangular';

Qx = 4;
Qy = 4;
Qc2 = 15;    % to see later

d = 1024;
flag_algo = algo_version;

% if lowRes
%     param_real_data.image_size_Nx = 2560; % 2560;
%     param_real_data.image_size_Ny = 2560; % 1536;
% else
%     param_real_data.image_size_Nx = 4096; % 2560;
%     param_real_data.image_size_Ny = 4096; % 1536;
% end
param_real_data.image_size_Nx = 2560;
param_real_data.image_size_Ny = 1536;
nChannels = length(ch);
% quite restrictive here), change l.70 accordingly
% klargestpercent = 20;

%% Config parameters
Nx = param_real_data.image_size_Nx;
Ny = param_real_data.image_size_Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 7; % number of neighbours for nufft
Ky = 7; % number of neighbours for nufft

if ~wterm
    [A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy * Ny ox * Nx], [Ny / 2 Nx / 2]);
end

%% Load data

for i = 1:nChannels
    fprintf('\nChannel number: %d\n', ch(i));
    DRfilename = [datadir, '/CYG_DR_cal_', num2str(realdatablocks), 'b_ind', ...
    num2str(subInd(1)), '_', num2str(subInd(end)), '_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(ch(i)), '.mat'];
%     if wterm
%         if lowRes
%             DRfilename = [datadir, '/ESO137_LOW_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '=', num2str(ch(i)), '.mat'];
%         else
%             DRfilename = [datadir, '/ESO137_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '=', num2str(ch(i)), '.mat'];
%         end
%         fprintf('Read dimensionality reduction file: %s\n', DRfilename)
%         tmp = load(DRfilename, 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'A', 'At', 'epsilon');
%     else
%         if lowRes
%             DRfilename = [datadir, '/ESO137_LOW_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(ch(i)), '.mat'];
%         else
%             DRfilename = [datadir, '/ESO137_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(ch(i)), '.mat'];
%         end
%         fprintf('Read dimensionality reduction file: %s\n', DRfilename)
%         tmp = load(DRfilename, 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'epsilon');
%     end
    tmp = load(DRfilename, 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'epsilon');
    H{i, 1} = tmp.H{1, 1};
    W{i, 1} = tmp.W{1, 1};
    yT{i, 1} = tmp.yT{1, 1};
    T{i, 1} = tmp.T{1, 1};
    aW{i, 1} = tmp.aW{1, 1};
    Wm{i, 1} = tmp.Wm{1, 1};
    epsilon{i, 1} = tmp.epsilon{1, 1};
    if wterm
        A = tmp.A;
        At = tmp.At;
    end
    % cheap H matrices
    if reduction_version == 2
        for j = 1:length(H{i, 1})
%         for j = 1:length(H)
%             H{i,1}{j} = H{i,1}{j}(Wm{i,1}{j}, :);
%             if adapt_eps_flag
%                 epsilon1{i,1}{j} = 0.9 * epsilon{i,1}{j};       % adjust epsilon for high frequency calibrated data
% %                 epsilon1{j} = 0.9 * epsilon{j};
%             end
            if usingPrecondition
                aW{i, 1}{j} = T{i, 1}{j};
%                 aW{j} = T{j};
            end
        end
    end
%     tmp = load(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_gam01_epsilon=', num2str(ch(i)), '.mat'], 'epsilon');
%     tmp = load(['/home/basphw/mjiang/Data/mjiang/real_data_dr/CYG_old_epsilon=', num2str(i), '.mat'], 'epsilon');
%     tmp = load(['/lustre/home/shared/sc004/dr_', num2str(realdatablocks), 'b_result_real_data/old/CYG_eps_cal_', num2str(realdatablocks), 'b_ind6_fouRed',...
%         num2str(reduction_version), '_th', num2str(fouRed_gamma),'=', num2str(ch(i)), '.mat'], 'epsilon');
%     tmp = load(['../../extract_real_data/CYG_eps_cal_', num2str(realdatablocks), 'b_ind6_fouRed',...
%         num2str(reduction_version), '_th', num2str(fouRed_gamma),'=', num2str(ch(i)), '.mat'], 'epsilon');
%     xsol(:,:,i) = tmp.xsol;
%     xsol = tmp.xsol;
    fprintf('\nDR file of channel number %d has been read\n', ch(i));

%     if wterm
%         DRfilename = [datadir, '/ESO137_NNLS_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '=', num2str(ch(i)), '.mat'];
%         fprintf("Read 'epsilon' from: %s\n", DRfilename)
%         tmp = load(DRfilename, 'epsilon');
%     else
%         DRfilename = [datadir, '/ESO137_NNLS_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(ch(i)), '.mat'];
%         fprintf("Read 'epsilon' from: %s\n", DRfilename)
%         tmp = load(DRfilename, 'epsilon');
%     end
%     epsilon{i,1} = tmp.epsilon{1,1};
end

clear tmp;

%% Compute full measurement operator spectral norm
if wterm
    if lowRes
        Anormfile = ['Anorm_LOW_dr_prec_ch', num2str(ch(1)), '_', num2str(ch(end)), '_ind', num2str(subInd(1)), '_', num2str(subInd(end)), '_', ...
        num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '.mat'];
    else
        Anormfile = ['Anorm_dr_prec_ch', num2str(ch(1)), '_', num2str(ch(end)), '_ind', num2str(subInd(1)), '_', num2str(subInd(end)), '_', ...
        num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '.mat'];
    end
else
    if lowRes
        Anormfile = ['Anorm_LOW_dr_prec_ch', num2str(ch(1)), '_', num2str(ch(end)), '_ind', num2str(subInd(1)), '_', num2str(subInd(end)), '_', ...
        num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.mat'];
    else
        Anormfile = ['Anorm_dr_prec_ch', num2str(ch(1)), '_', num2str(ch(end)), '_ind', num2str(subInd(1)), '_', num2str(subInd(end)), '_', ...
        num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.mat'];
    end
end
if isfile(Anormfile)
    compute_Anorm = false;
else
    compute_Anorm = true;
end

if compute_Anorm
    fprintf('\nCompute the operator norm: \n');
    if reduction_version == 1
        F = afclean(@(x) HS_fouRed_forward_operator(x, A, At, H, W, T, Wm, aW));
        Ft = afclean(@(y) HS_fouRed_adjoint_operator(y, A, At, H, W, T, Wm, aW));
    elseif reduction_version == 2
        F = afclean(@(x) HS_operatorGtPhi(x, H, W, A, T, aW));
        Ft = afclean(@(y) HS_operatorGtPhi_t(y, H, W, At, T, aW));
    end
    Anorm = pow_method_op(F, Ft, [Ny Nx nChannels]);
    save(Anormfile, '-v7.3', 'Anorm');
else
    fprintf('\nLoad the operator norm file: %s\n', Anormfile);
    load(Anormfile);
end
% delete(gcp('nocreate'));
fprintf('\nThe operator norm: %e\n', Anorm);

clear F Ft;

%% Sparsity operator definition
nlevel = 4; % wavelet level
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8'}; % 'self' % wavelet basis to be used
L = [2 * (1:8)'; 0]; % length of the filters (0 corresponding to the 'self' basis)

if flag_algo < 2

    [Psi1, Psit1] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    P = length(Psi1);

    for k = 1:P
        f = '@(x_wave) HS_forward_sparsity(x_wave,Psi1{';
        f = sprintf('%s%i},Ny,Nx);', f, k);
        Psi{k} = eval(f);

        b(k) = size(Psit1{k}(zeros(Ny, Nx, 1)), 1);

        ft = ['@(x) HS_adjoint_sparsity(x,Psit1{' num2str(k) '},b(' num2str(k) '));'];
        Psit{k} = eval(ft);
    end

    %% Full sparsity operator
    [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    bb = size(Psitw(zeros(Ny, Nx, 1)), 1);

    Psi_full = @(x_wave) HS_forward_sparsity(x_wave, Psiw, Ny, Nx);
    Psit_full = @(x) HS_adjoint_sparsity(x, Psitw, bb);

elseif flag_algo == 2
    [Psi, Psit] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
end

%% L21 + Nuclear (facet-based version)
if algo_version == 1
    disp('Split L21 + Nuclear + wavelets');
    %% HSI parameter structure sent to the  HSI algorithm
    param_HSI.verbose = 2; % print log or not
    param_HSI.nu0 = 1; % bound on the norm of the Identity operator
    param_HSI.nu1 = 1; % bound on the norm of the operator Psi
    param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
    param_HSI.gamma0 = gamma0;
    param_HSI.gamma = gamma;  % convergence parameter L1 (soft th parameter)
    param_HSI.rel_obj = 1e-10; % stopping criterion
    param_HSI.max_iter = 10; % max number of iterations

    param_HSI.use_adapt_eps = adapt_eps_flag; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
    param_HSI.adapt_eps_start = 300; % minimum num of iter before stating adjustment
    param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
    param_HSI.adapt_eps_tol_out = 1.01; % tolerance outside the l2 ball
    param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
    param_HSI.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
    param_HSI.adapt_eps_change_percentage = 0.5 * (sqrt(5) - 1); % the weight of the update w.r.t the l2 norm of the residual data
%     param_HSI.l2_upper_bound = epsilon;

    param_HSI.reweight_alpha = (0.8)^10; % 1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
    param_HSI.reweight_alpha_ff = 0.8;
    param_HSI.total_reweights = 50; % -1 if you don't want reweighting
    param_HSI.reweight_abs_of_max = Inf; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised

    param_HSI.use_reweight_steps = 0; % reweighting by fixed steps
    param_HSI.reweight_step_size = 300; % reweighting step size
    param_HSI.reweight_steps = [5000:param_HSI.reweight_step_size:10000];
    param_HSI.step_flag = 1;

    param_HSI.use_reweight_eps = 1; % reweighting w.r.t the relative change of the solution
    param_HSI.reweight_max_reweight_itr = param_HSI.max_iter - param_HSI.reweight_step_size;
    param_HSI.reweight_rel_obj = 5e-4; % criterion for performing reweighting
    param_HSI.reweight_min_steps_rel_obj = 300; % min num of iter between reweights

    param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid
    param_HSI.precondition = usingPrecondition;
    param_HSI.rw_tol = 5000;
    param_HSI.ind = 1:16;

%     param_HSI.initsol = xsol;

    reweight_step_count = 0;
    initfilename = init_file_name;

    % spectral tesselation (non-overlapping)
    rg_c = domain_decomposition(Qc2, nChannels); % to be modified % only for data reduction
    cell_c_chunks = cell(Qc2, 1);
    y_spmd = cell(Qc2, 1);
    epsilon_spmd = cell(Qc2, 1);
    l2_upper_bound_spmd = cell(Qc2, 1);
    aW_spmd = cell(Qc2, 1);
    Wm_spmd = cell(Qc2, 1);
    T_spmd = cell(Qc2, 1);
    H_spmd = cell(Qc2, 1);
    W_spmd = cell(Qc2, 1);

    for i = 1:Qc2
        cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
%         cell_c_chunks{i} = [rg_c(i, 1):rg_c(i, 2), nChannels-rg_c(i, 2)+1:nChannels-rg_c(i, 1)+1];   % only for data reduction
        y_spmd{i} = yT(cell_c_chunks{i});
        if adapt_eps_flag
            epsilon_spmd{i} = epsilon(cell_c_chunks{i});
        else
            epsilon_spmd{i} = epsilon(cell_c_chunks{i});
        end
        l2_upper_bound_spmd{i} = epsilon(cell_c_chunks{i});
        aW_spmd{i} = aW(cell_c_chunks{i});
        Wm_spmd{i} = Wm(cell_c_chunks{i});
        T_spmd{i} = T(cell_c_chunks{i});
        H_spmd{i} = H(cell_c_chunks{i});
        W_spmd{i} = W(cell_c_chunks{i});
    end
    param_HSI.l2_upper_bound = l2_upper_bound_spmd;
    clear yT epsilon aW Wm T epsilon;

    if  rw >= 0
        load(['results/result_HyperSARA_spmd4_cst_weighted_rd_' num2str(param_HSI.gamma) '_' num2str(rw) '.mat']);
        param_HSI.init_xsol = param.init_xsol;
        param_HSI.init_g = param.init_g;
        param_HSI.init_v0 = param.init_v0;
        param_HSI.init_v1 = param.init_v1;
        param_HSI.init_weights0 = param.init_weights0;
        param_HSI.init_weights1 = param.init_weights1;
        param_HSI.init_v2 = param.init_v2;
        param_HSI.init_proj = param.init_proj;
        param_HSI.init_t_block = param.init_t_block;
        param_HSI.init_t_start = param.init_t_start + 1;
        param_HSI.reweight_alpha = param.reweight_alpha;
        param_HSI.init_reweight_step_count = param.init_reweight_step_count;
        param_HSI.init_reweight_last_iter_step = param.init_reweight_last_iter_step;
        param_HSI.reweight_steps = (param_HSI.init_t_start + param_HSI.reweight_step_size:param_HSI.reweight_step_size:param_HSI.max_iter + (2 * param_HSI.reweight_step_size));
        param_HSI.step_flag = 0;
    end

    % solvers
    mkdir('results/');
    [xsol, param_HSI, t, rel_fval, nuclear, l21, norm_res_out, res, end_iter] = ...
        facetHyperSARA_DR_precond(y_spmd, epsilon_spmd, A, At, H_spmd, W_spmd, aW_spmd, T_spmd, Wm_spmd, param_HSI, Qx, Qy, Qc2, ...
        wlt_basis, L, nlevel, cell_c_chunks, nChannels, d, window_type, initfilename, reduction_version, realdatablocks, fouRed_gamma, typeStr);

    save(['results/results_facethyperSARA_fouRed_ch', num2str(ch(1)), '_', num2str(ch(end)), '_', num2str(algo_version), '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), ...
        '_Qc=', num2str(Qc2), '_gamma=', num2str(gamma), '_gamma0=', num2str(gamma0), '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_adpteps', num2str(adapt_eps_flag), '.mat'], '-v7.3', ...
        'xsol', 'param_HSI', 't', 'rel_fval', 'nuclear', 'l21', 'norm_res_out', 'res', 'end_iter');

elseif algo_version == 2
    disp('Split L1 + wavelets');
    % %% PDFB parameter structure sent to the algorithm
    param_pdfb.verbose = 2; % print log or not
    param_pdfb.nu1 = 1; % bound on the norm of the operator Psi
    param_pdfb.nu2 = Anorm; % bound on the norm of the operator A*G
    param_pdfb.gamma = gamma; % convergence parameter L1 (soft th parameter)
    param_pdfb.tau = 0.49; % forward descent step size
    param_pdfb.rel_obj = 1e-10; % stopping criterion
    param_pdfb.max_iter = 10000; % max number of iterations
    param_pdfb.lambda0 = 1; % relaxation step for primal update
    param_pdfb.lambda1 = 1; % relaxation step for L1 dual update
    param_pdfb.lambda2 = 1; % relaxation step for L2 dual update
    param_pdfb.omega1 = 1;
    param_pdfb.omega2 = 1;
    param_pdfb.sol_steps = [inf]; % saves images at the given iterations

    param_pdfb.use_proj_elipse_fb = logical(usingPrecondition);
    param_pdfb.elipse_proj_max_iter = 20;
    param_pdfb.elipse_proj_min_iter = 1;
    param_pdfb.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid

    param_pdfb.use_adapt_eps = adapt_eps_flag; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
    param_pdfb.adapt_eps_start = 300; % minimum num of iter before stating adjustment
    param_pdfb.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
    param_pdfb.adapt_eps_tol_out = 1.001; % tolerance outside the l2 ball
    param_pdfb.adapt_eps_steps = 100; % min num of iter between consecutive updates
    param_pdfb.adapt_eps_rel_obj = 1e-3; % bound on the relative change of the solution
    param_pdfb.adapt_eps_change_percentage = 0.5 * (sqrt(5) - 1); % the weight of the update w.r.t the l2 norm of the residual data
    param_pdfb.l2_upper_bound = epsilon{1};

    param_pdfb.reweight_alpha = 0.8^10; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
    param_pdfb.reweight_alpha_ff = rw_alpha;
    param_pdfb.total_reweights = rw_tot; % -1 if you don't want reweighting
    param_pdfb.reweight_abs_of_max = Inf; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised

    param_pdfb.use_reweight_steps = 0; % reweighting by fixed steps
    param_pdfb.reweight_step_size = 300; % reweighting step size
    param_pdfb.reweight_steps = [5000:param_pdfb.reweight_step_size:10000];
    param_pdfb.step_flag = 1;

    param_pdfb.use_reweight_eps = 1; % reweighting w.r.t the relative change of the solution
    param_pdfb.reweight_max_reweight_itr = param_pdfb.max_iter - param_pdfb.reweight_step_size;
    param_pdfb.reweight_rel_obj = 1e-4; % criterion for performing reweighting
    param_pdfb.reweight_min_steps_rel_obj = 100; % min num of iter between reweights

    param_pdfb.ppd_max_iter = 300;
    param_pdfb.rw_tol = 5000;

%     param_pdfb.initsol = xsol;
    param_pdfb.chInd = ch(1);
    param_pdfb.init_file_name = init_file_name;

    % solvers
    mkdir('results/');
    [xsol, ~, epsilon, t, rel_fval, norm2, res, end_iter] = ...
        pdfb_DR_precond(yT{1}, epsilon{1}, A, At, H{1}, W{1}, aW{1}, T{1}, Wm{1}, ...
        Psi, Psit, param_pdfb, reduction_version, realdatablocks, fouRed_gamma, typeStr, wterm, levelG, levelC, init_file_name);

    save(['results/results_fouRed_ch', num2str(ch(1)), '_', num2str(ch(end)), '_ind', num2str(subInd(1)), '_', num2str(subInd(end)), ...
        '_', num2str(algo_version), '_gamma=', num2str(gamma), '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.mat'], '-v7.3', ...
        'xsol', 'epsilon', 't', 'rel_fval', 'norm2', 'res', 'end_iter');

end

end
