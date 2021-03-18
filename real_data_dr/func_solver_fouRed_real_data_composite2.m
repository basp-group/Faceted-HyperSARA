function func_solver_fouRed_real_data_composite2(datadir, name, Qx, Qy, Qc2, gamma0, gamma1, ch, subInd, reduction_version, algo_version, realdatablocks, fouRed_gamma, fouRed_type, adapt_eps_flag, jobpath, flag_primal, flag_homotopy)
    % Global imaging solver for DR real data 
    %   version "composite", parcluster initialised outside the main solver
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
    % > reduction_version: reduction version, 
    %       1: old one (not used any more)
    %       2: current one
    % > algo_version: version of algorithm
    %       1. hyperspectral version
    % > realdatablocks: number of data blocks
    % > fouRed_gamma: level of reduction
    % > fouRed_type: type of reduction, supplementary information for
    % fouRed_gamma
    %       1: remove smallest singular values based on "fouRed_gamma" percentage 
    %       2: remove smallest singular values based on "fouRed_gamma"-sigma 
    % > adapt_eps_flag: flag of adaptive epsilon strategy
    % > jobpath: path of job (for runnings on cluster)
    
    if fouRed_type == 1
        typeStr = 'perc';
    elseif fouRed_type == 2
        typeStr = 'th';
    end
    
    addpath ../fouRed
    addpath ../lib/
    addpath ../lib/operators/
    addpath ../lib/measurement-operator/nufft/
    addpath ../lib/utils/
    addpath ../lib/faceted-wavelet-transform/src
    addpath ../src_mnras/
    addpath ../src_mnras/spmd/
    addpath ../src_mnras/spmd/dr/
    addpath ../src_mnras/spmd/weighted/
    
    fprintf("Nuclear norm parameter=%e\n", gamma0);
    fprintf("Regularization parameter=%e\n", gamma1);
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
    
    window_type = 'triangular';
    
    d = 512;
    flag_algo = algo_version;
    
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
    tau_ch = 0.5;
    nChannels1 = tau_ch*nChannels;
    nChannels2 = (1-tau_ch)*nChannels;
    rg_c1 = split_range(Qc2, nChannels1); % only for data reduction
    rg_c2 = split_range(Qc2, nChannels2); % only for data reduction
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
    cirrus_cluster = parcluster('local'); % 'local' 'SlurmProfile1'
    cirrus_cluster.JobStorageLocation= jobpath;
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
    filter_length = [2*(1:8)'; 0]; % length of the filters (0 corresponding to the 'self' basis)

    %! -- TO BE CHECKED
    % compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
    % SARA space) involved in the reweighting scheme
    [sig, sig_bar, max_psf, ~, ~, ~] = ...
    compute_reweighting_lower_bound_dr(yTp, Wp, Tp, Hp, Ap, Atp, Ny, Nx, ...
    nChannels, wlt_basis, filter_length, nlevel, Q, cell_c_chunks);
    %! --

    %% L21 + Nuclear (facet-based version)
    if flag_algo == 1
        disp('Split L21 + Nuclear + wavelets')    
        %% HSI parameter structure sent to the  HSI algorithm
        param_HSI.verbose = 2; % print log or not
        param_HSI.nu0 = 1; % bound on the norm of the Identity operator
        param_HSI.nu1 = 1; % bound on the norm of the operator Psi
        param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
        param_HSI.gamma0 = gamma0;
        param_HSI.gamma1 = gamma1;  %convergence parameter L1 (soft th parameter)
        param_HSI.rel_val = 1e-10; % stopping criterion
        param_HSI.max_iter = 100000; % max number of iterations
    
        param_HSI.use_adapt_eps = adapt_eps_flag; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
        param_HSI.adapt_eps_start = 300; % minimum num of iter before stating adjustment
        param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
        param_HSI.adapt_eps_tol_out = 1.01; % tolerance outside the l2 ball
        param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
        param_HSI.adapt_eps_rel_var = 5e-4; % bound on the relative change of the solution
        param_HSI.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data
        
        %! -- TO BE CHECKED: see where reweight_alpha needs to start from        
        if flag_homotopy
            param_HSI.reweight_alpha = 10;
            param_HSI.reweight_alpha_ff = (1/param_HSI.reweight_alpha)^(1/6); % reach the floor level in 6 reweights (see if a different number would be appropriate)
            % param_HSI.reweight_alpha = (0.8)^10; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
            % param_HSI.reweight_alpha_ff = 0.8;
        else
            param_HSI.reweight_alpha = 1;
            param_HSI.reweight_alpha_ff = 1;
        end
        %! --

        param_HSI.total_reweights = 50; % -1 if you don't want reweighting
        param_HSI.reweight_abs_of_max = Inf; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised
        param_HSI.sig = sig; % estimate of the noise level in SARA space
        param_HSI.sig_bar = sig_bar; % estimate of the noise level in "SVD" space
        param_HSI.max_psf = max_psf;

        param_HSI.use_reweight_steps = 0; % reweighting by fixed steps
        param_HSI.reweight_step_size = 300; % reweighting step size
        param_HSI.reweight_steps = [5000: param_HSI.reweight_step_size :10000];
        param_HSI.step_flag = 1;
    
        param_HSI.use_reweight_eps = 1; % reweighting w.r.t the relative change of the solution
        param_HSI.reweight_max_reweight_itr = param_HSI.max_iter - param_HSI.reweight_step_size;
        param_HSI.reweight_rel_val = 5e-4; % criterion for performing reweighting
        param_HSI.reweight_min_steps_rel_obj = 300; % min num of iter between reweights
    
        param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
        param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
        param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid
        param_HSI.precondition = usingPrecondition;
        param_HSI.ind = 1:16;
        
        param_HSI.ppd_min_iter = 500;
        param_HSI.ppd_max_iter = 500;
        
        param_HSI.initsol = xsol;
        
        reweight_step_count = -1;
        initfilename = ['./results/', name, '_dr_co_w_real_' ...
                    num2str(param_HSI.ind(1)), '_', num2str(param_HSI.ind(end)), '_', num2str(param_HSI.gamma1), '_' num2str(param_HSI.gamma0), '_adpteps', num2str(param_HSI.use_adapt_eps), '_' num2str(reweight_step_count) '_' ... 
                    '_primal=' num2str(flag_primal), '_homotopy=', num2str(flag_homotopy) '.mat'];
            
        % spectral tesselation (non-overlapping)
        epsilon_spmd = cell(Qc2, 1);
           
        for i = 1:Qc2
            if adapt_eps_flag
                epsilon_spmd{i} = epsilon(cell_c_chunks{i});
            else
                epsilon_spmd{i} = epsilon(cell_c_chunks{i});
            end
            l2_upper_bound{i} = epsilon(cell_c_chunks{i});
        end
        
        param_HSI.l2_upper_bound = l2_upper_bound;
        
        clear epsilon epsilon1
        
        % solvers
        mkdir('results/')
    
        [xsol,param_HSI,t,rel_fval,nuclear,l21,norm_res_out,end_iter] = ...
        facetHyperSARA_DR_precond_v22(yTp, epsilon_spmd, Ap, Atp, Hp, Wp, aWp, Tp, Wmp, param_HSI, ...
        Qx, Qy, Qc2, wlt_basis, filter_length, nlevel, cell_c_chunks, nChannels, d, window_type, initfilename, name, ...
        reduction_version, realdatablocks, fouRed_gamma, typeStr, Ny, Nx, flag_primal, flag_homotopy);
        
        save(['results/results_', name, '_', num2str(ch(1)), '_', num2str(ch(end)), '_', num2str(algo_version), '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), ...
            '_Qc=', num2str(Qc2), '_gamma=', num2str(gamma1), '_gamma0=', num2str(gamma0), '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_adpteps', num2str(adapt_eps_flag),'_primal=' num2str(flag_primal), '_homotopy=', num2str(flag_homotopy),'.mat'], '-v7.3', ...
            'xsol', 'param_HSI', 't', 'rel_fval', 'nuclear', 'l21', 'norm_res_out', 'end_iter');
    end
