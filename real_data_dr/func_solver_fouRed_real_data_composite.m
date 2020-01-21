function func_solver_fouRed_real_data_composite(datadir, gamma0, gamma, ch, subInd, reduction_version, algo_version, realdatablocks, fouRed_gamma, fouRed_type, adapt_eps_flag)

if fouRed_type == 1
    typeStr = 'perc';
elseif fouRed_type == 2
    typeStr = 'th';
end

diaryFname = ['diary_ch', num2str(ch(1)), '_', num2str(ch(end)), '_ind', num2str(subInd(1)), '_', num2str(subInd(end)), ...
    '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_algo', num2str(algo_version), '_', typeStr, num2str(fouRed_gamma), '.txt'];
    
if exist(diaryFname, 'file')
    delete(diaryFname)
end
diary(diaryFname)

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

% compute_Anorm = true;
usingPrecondition = true;
rw = -1;
window_type = 'triangular';

Qx = 5;
Qy = 3;
Qc2 = 6;    % to see later

d = 512;
flag_algo = algo_version;
% parallel_version = 'spmd4_cst';
% bool_weights = true; % for the spmd4_new version (50% overlap version)

param_real_data.image_size_Nx = 2560; % 2560;
param_real_data.image_size_Ny = 1536; % 1536;
nChannels = length(ch); % total number of "virtual" channels (i.e., after
% concatenation) for the real dataset considered
% nBlocks = realdatablocks;        % number of data blocks (needs to be known beforehand,
% quite restrictive here), change l.70 accordingly
% klargestpercent = 20;

%% Config parameters
Nx = param_real_data.image_size_Nx;
Ny = param_real_data.image_size_Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 7; % number of neighbours for nufft
Ky = 7; % number of neighbours for nufft

[A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);

% spectral tesselation (non-overlapping)
tau_ch = 0.6;
nChannels1 = tau_ch*nChannels;
nChannels2 = (1-tau_ch)*nChannels;
rg_c1 = domain_decomposition(Qc2, nChannels1); % to be modified % only for data reduction
rg_c2 = domain_decomposition(Qc2, nChannels2); % to be modified % only for data reduction
cell_c_chunks = cell(Qc2, 1);
nchannel_per_worker = zeros(Qc2, 1);

for i = 1:Qc2
%     cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
    cell_c_chunks{i} = [rg_c1(i, 1):rg_c1(i, 2), nChannels-rg_c2(i, 2)+1:nChannels-rg_c2(i, 1)+1];   % only for data reduction
    nchannel_per_worker(i) = numel(cell_c_chunks{i});
end
    
%% Load data
delete(gcp('nocreate'));
Q = Qx*Qy; 
numworkers = Q + Qc2;
cirrus_cluster = parcluster('cirrus R2019a');
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
        for i = 1: ch_len
            ch_ind = chunk(i);
            fprintf('\nChannel number: %d\n', ch(ch_ind))
        %     DRfilename = ['/lustre/home/shared/sc004/dr_', num2str(realdatablocks), 'b_result_real_data/CYG_DR_cal_', num2str(realdatablocks), 'b_ind',...
        %     num2str(subInd(1)), '_', num2str(subInd(end)), '_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma),'=', num2str(ch(i)), '.mat'];
            DRfilename = [datadir, '/CYG_DR_cal_', num2str(realdatablocks), 'b_ind', num2str(subInd(1)), '_', num2str(subInd(end)), '_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma),'=', num2str(ch(ch_ind)), '.mat'];
            fprintf('Read dimensionality reduction file: %s\n', DRfilename)
            tmp = load(DRfilename, 'H', 'W', 'yT', 'T', 'aW', 'Wm');
            Hp{i,1} = tmp.H{1,1};
            Wp{i,1} = tmp.W{1,1};
            yTp{i,1} = tmp.yT{1,1};
            Tp{i,1} = tmp.T{1,1};
            aWp{i,1} = tmp.aW{1,1};
            Wmp{i,1} = tmp.Wm{1,1};
            if reduction_version == 2
                for j = 1:length(Hp{i,1})
                    if usingPrecondition
                        aWp{i,1}{j} = Tp{i,1}{j};
                    end
                end
            end
            fprintf('\nDR file of channel number %d has been read\n', ch(ch_ind))
        end
    end
end

clear tmp
 
for i = 1:nChannels
    DRfilename = [datadir, '/CYG_DR_cal_', num2str(realdatablocks), 'b_ind', num2str(subInd(1)), '_', num2str(subInd(end)), '_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma),'=', num2str(ch(i)), '.mat'];
    fprintf("Read 'epsilon', 'xsol' from: %s\n", DRfilename)
    tmp = load(DRfilename, 'epsilon', 'xsol');
    epsilon{i,1} = tmp.epsilon{1,1};
        % cheap H matrices
        if reduction_version == 2
            for j = 1:length(epsilon{i,1})
                if adapt_eps_flag
                    epsilon1{i,1}{j} = 0.9 * epsilon{i,1}{j};       % adjust epsilon for high frequency calibrated data
                end
            end
        end
    xsol(:,:,i) = tmp.xsol;
end

%% Compute full measurement operator spectral norm
Anormfile = ['Anorm_dr_prec_ch', num2str(ch(1)), '_', num2str(ch(end)), '_ind', num2str(subInd(1)), '_', num2str(subInd(end)), '_',...
    num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.mat'];
if isfile(Anormfile)
    compute_Anorm = false;
else
    compute_Anorm = true;
end

if compute_Anorm
    fprintf('\nCompute the operator norm: \n')
%     spmd
%         if labindex > Q
%             if reduction_version == 2
%                 Fp = afclean( @(x) HS_operatorGtPhi(x, Hp, Wp, Ap, Tp, aWp));
%                 Ftp = afclean( @(y) HS_operatorGtPhi_t(y, Hp, Wp, Atp, Tp, aWp));
%             end
%         end
%     end
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

%% L21 + Nuclear (facet-based version)
if algo_version == 1
    disp('Split L21 + Nuclear + wavelets')    
    %% HSI parameter structure sent to the  HSI algorithm
    param_HSI.verbose = 2; % print log or not
    param_HSI.nu0 = 1; % bound on the norm of the Identity operator
    param_HSI.nu1 = 1; % bound on the norm of the operator Psi
    param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
    param_HSI.gamma0 = gamma0;
    param_HSI.gamma = gamma;  %convergence parameter L1 (soft th parameter)
    param_HSI.rel_obj = 1e-10; % stopping criterion
    param_HSI.max_iter = 10000; % max number of iterations

    param_HSI.use_adapt_eps = adapt_eps_flag; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
    param_HSI.adapt_eps_start = 300; % minimum num of iter before stating adjustment
    param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
    param_HSI.adapt_eps_tol_out = 1.01; % tolerance outside the l2 ball
    param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
    param_HSI.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
    param_HSI.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data

    param_HSI.reweight_alpha = (0.8)^10; %1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
    param_HSI.reweight_alpha_ff = 0.8;
    param_HSI.total_reweights = 50; % -1 if you don't want reweighting
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
    param_HSI.rw_tol = 5000;
    param_HSI.ind = 1:16;
    
    param_HSI.initsol = xsol;
    
    reweight_step_count = 0;
    initfilename = ['./results/facetHyperSARA_dr_co_w_real_' ...
                num2str(param_HSI.ind(1)), '_', num2str(param_HSI.ind(end)), '_' num2str(param_HSI.gamma) '_' num2str(reweight_step_count) '.mat'];
        
    % spectral tesselation (non-overlapping)
    epsilon_spmd = cell(Qc2, 1);
       
    for i = 1:Qc2
        if adapt_eps_flag
            epsilon_spmd{i} = epsilon1(cell_c_chunks{i});
        else
            epsilon_spmd{i} = epsilon(cell_c_chunks{i});
        end
        l2_upper_bound{i} = epsilon(cell_c_chunks{i});
    end
    
    param_HSI.l2_upper_bound = l2_upper_bound;
    
    clear epsilon epsilon1
    
    if  rw >= 0 
        load(['results/result_HyperSARA_spmd4_cst_weighted_rd_' num2str(param_HSI.gamma) '_' num2str(rw) '.mat']);
        
        if adapt_eps_flag
            epsilon_spmd = epsilon1;
        else
            epsilon_spmd = epsilon;
        end
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
    
    [xsol,param_HSI,t,rel_fval,nuclear,l21,norm_res_out,end_iter] = ...
    facetHyperSARA_DR_precond_v2(yTp, epsilon_spmd, Ap, Atp, Hp, Wp, aWp, Tp, Wmp, param_HSI, ...
    Qx, Qy, Qc2, wlt_basis, L, nlevel, cell_c_chunks, nChannels, d, window_type, initfilename, ...
    reduction_version, realdatablocks, fouRed_gamma, typeStr, Ny, Nx);
    
    save(['results/results_facethyperSARA_fouRed_ch', num2str(ch(1)), '_', num2str(ch(end)), '_', num2str(algo_version), '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), ...
        '_Qc=', num2str(Qc2), '_gamma=', num2str(gamma), '_gamma0=', num2str(gamma0), '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_adpteps', num2str(adapt_eps_flag),'.mat'], '-v7.3', ...
        'xsol', 'param_HSI', 't', 'rel_fval', 'nuclear', 'l21', 'norm_res_out', 'end_iter');
diary off
end
