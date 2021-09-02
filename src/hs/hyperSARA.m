function xsol = hyperSARA(y, epsilon, A, At, pU, G, W, param, K, ... 
    wavelet, nlevel, spectral_chunk, nChannels, M, N, oy, ox, ...
    name_warmstart, name_checkpoint, flagDR, Sigma, varargin)

% TODO: to distribute before entering the solver: y, G, pU, W, X0
% TODO: try to replace A, At, pU, G, W by a functor (if possible)

%HyperSARA
%
% ...
%
%-------------------------------------------------------------------------%
%%
% Input: 
%
% > y           blocks of visibilities {L}{nblocks_l}
% > epsilon     l2-ball norms {L}{nblocks_l}
% > A           measurement operator
% > At          adjoint of the measurement operator
% > pU          preconditioning matrices {L}{nblocks_l}
% > G           gridding matrices {L}{nblocks_l}
% > W           masks for selection of the blocks of visibilities
% > param       algorithm parameters (struct)
% > X0          ground truth image cube
% > K           number of data workers
% > wavelet     list of wavelet basis (for the SARA prior)
% > nlevel      number of wavelet decomposition levels
% > spectral_chunk    list of channels handled by each data process
% > nChannels           total number of channels
% > name_warmstart name_checkpoint of the file to restart from
% > name_checkpoint        lambda function defining the name_checkpoint of the backup file 
% > flag_homotopy flag to activate homotopy scheme in the reweighting scheme
% > varargin     initial value for the primal variable
%
%   general
%   > .verbose  print log or not
%
%   convergence
%   > .nu0 = 1
%   > .nu1      upper bound on the norm of the operator Psi
%   > .nu2      upper bound on the norm of the measurement operator
%   > .gamma0   regularization parameter (nuclear norm)
%   > .gamma    regularization parameter (l21 norm)
%
%   pdfb
%   > .pdfb_min_iter               minimum number of iterations
%   > .pdfb_max_iter               maximum number of iterations
%   > .pdfb_rel_var                relative variation tolerance
%   > .pdfb_fidelity_tolerance     tolerance to check data constraints are satisfied
%
%   reweighting       
%   > .reweighting_max_iter  (30)    maximum number of reweighting steps %! (weights updated reweighting_max_iter - 1 times)
%   > .reweighting_min_iter          minimum number of reweighting steps (to reach "noise" level)
%   > .reweighting_rel_var  (1e-4)   relative variation
%   > .reweighting_alpha             starting reweighting parameter (> 1)
%   > .reweighting_alpha_ff  (0.9)   multiplicative parameter update (< 1) 
%   > .reweighting_sig               noise level (wavelet space)
%   > .reweighting_sig_bar           noise level (singular value space)
%
%   projection onto ellipsoid (preconditioning)
%   > .elipse_proj_max_iter (20)     max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
%   > .elipse_proj_min_iter  (1)     min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
%   > .elipse_proj_eps    (1e-8)     stopping criterion
%
%   adaptive epsilon
%   > .use_adapt_eps                 flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
%   > .adapt_eps_start (200)         minimum num of iter before stating adjustment
%   > .adapt_eps_tol_in (0.99)       tolerance inside the l2 ball
%   > .adapt_eps_tol_out (1.001)     tolerance outside the l2 ball
%   > .adapt_eps_steps (100)         min num of iter between consecutive updates
%   > .adapt_eps_rel_var (5e-4)      bound on the relative change of the solution
%   > .adapt_eps_change_percentage  (0.5*(sqrt(5)-1)) weight of the update w.r.t the l2 norm of the residual data
%
% > X0          ground truth wideband image [M*N, L]
% > K           number of Matlab data fidelity processes [1]
% > wavelet     wavelet doctionaries considered (should contain 'self' by
%               default in last position)
% > nlevel      decomposition depth [1]
% > spectral_chunk    indices of the bands handled by each data node {K, 1}
% > nChannels           total number of spectral channels [1]
% > name_warmstart  name_checkpoint of a valid .mat file for initialization (for warm-restart)
% > name_checkpoint        lambda function defining the name_checkpoint of the backup file 
% > flag_homotopy flag to activate homotopy scheme in the reweighting scheme
% > varargin     initial value for the primal variable
%
%
% Output:
%
% < xsol        reconstructed wideband image [M*N, L]
% < v0          dual variables associated with the nuclear norms {Q}
% < v1          dual variables associated with the l21 norms {Q}
% < v2          dual variables associated with the data fidelity terms {K}
% < weights0    weights associated with the nuclear norm prior {Q}
% < weights1    weights associated with the nuclear norm prior {Q}
% < proj        projected 
% < t_block     index of the last iteration where the weigths have been
%               updated
% < reweighting_alpha  last value of the reweigthing parameter [1]
% < epsilon         updated value of th l2-ball radii {...}
% < t               index of the last iteration step [1]
% < rel_val        relative variation
% < nuclear         value of the faceted nuclear norm
% < l21             value of the l21 regularization term
% < norm_res_out    norm of the reidual image 
% < res             residual image [M, N]
% < end_iter        last iteration
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin, A. Abdulaziz, M. Jiang
% [../../2019]
%-------------------------------------------------------------------------%
%%
% Note:
% Code based on the HyperSARA code developed by A. Abdulaziz, available at 
% https://basp-group.github.io/Hyper-SARA/
%-------------------------------------------------------------------------%
%%
%SPMD version: use spmd for all the priors, deal with the data fidelity
% term in a single place.

% This function solves:
%
% min || X ||_* + lambda * ||Psit(X)||_2,1   s.t.  || Y - A(X) ||_2 <= epsilon and x>=0
%
% Author: Abdullah Abdulaziz
% Modified by: P.-A. Thouvenin.

%% REMARKS
% 1. Do I need to have the number of nodes Q, I available as parallel
% constants? (automatically broadcasted otherwise each time they are used?)
% 2. Try to find a compromise between Li and Nk, depending on the size of
% the problem at hand (assuming the data size is not the main bottleneck),
% and depending on the computational complexity of the task to be driven by
% each worker -> see if there is any acceleration here...
% 3. Version w/o data blocks (not in use here...)
%%

% initialize monitoring variables (display active)
norm_epsilon_check = Inf;
norm_residual_check = 0;

% size of the oversampled Fourier space (vectorized)
No = M*oy*N*ox;

% instantiate Psi, Psit
spmd
    [Psi_, Psit_] = op_sp_wlt_basis(wavelet, nlevel, M, N);
end

% Initializations
init_flag = isfile(name_warmstart);
if init_flag
    init_m = matfile(name_warmstart);
    fprintf('Resume from file %s\n\n', name_warmstart)
end

%! -- TO BE CHECKED (primal initialization)
if init_flag
    xsol = init_m.xsol;
    pdfb_rel_var_low = param.pdfb_rel_var_low;
    param = init_m.param;

    if ~isfield(param,'pdfb_rel_var_low')
        param.pdfb_rel_var_low = pdfb_rel_var_low;
    end

    % check flag homotopy strategy
    if ~isfield(param, 'flag_homotopy')
        flag_homotopy = false;
    else
        flag_homotopy = param.flag_homotopy;
    end

    epsilon = Composite();
    for k = 1:K
        epsilon{k} = init_m.epsilon(spectral_chunk{k}, 1);
    end

    if numel(varargin) > 1
        flag_synth_data = true;
        X0 = varargin{2};
    else
        flag_synth_data = false;
    end
    fprintf('xsol, param and epsilon uploaded \n\n')

else
    if ~isempty(varargin)
        if ~isempty(varargin{1})
            xsol = varargin{1};
        else
            xsol = zeros(M, N, nChannels);
        end
        
        if numel(varargin) > 1
            flag_synth_data = true;
            X0 = varargin{2};
        else
            flag_synth_data = false;
        end
    else
        xsol = zeros(M, N, nChannels);
    end
    fprintf('xsol initialized \n\n')
end
%! assumes backup file exactly saved at the time a reweighting step occurred
xlast_reweight = xsol;
%! --

if init_flag
    g = init_m.g;
    fprintf('g uploaded \n\n')
else
    g = zeros(size(xsol));
    fprintf('g initialized \n\n')
end  

%! -- TO BE CHECKED
% Reweighting parameters
sig_bar_ = param.reweighting_sig_bar;
sig_ = Composite();
for k = 1:K
    sig_{k} = param.reweighting_sig;
end
reweighting_alpha = param.reweighting_alpha;
reweighting_alphap = Composite();
for k = 1:K
    reweighting_alphap{k} = reweighting_alpha;
end
% reweighting_alpha_ffp = parallel.pool.Constant(param.reweighting_alpha_ff);

if isfield(param,'init_reweight_step_count')
    reweight_step_count = param.init_reweight_step_count;
    fprintf('reweight_step_count uploaded\n\n')
else
    reweight_step_count = 0;
    fprintf('reweight_step_count initialized \n\n')
end

if isfield(param,'init_reweight_last_iter_step')
    reweight_last_step_iter = param.init_reweight_last_iter_step;
    fprintf('reweight_last_iter_step uploaded \n\n')
else
    param.init_reweight_last_iter_step = 0;
    reweight_last_step_iter = 0;
    fprintf('reweight_last_iter_step initialized \n\n')
end
%! --

% Primal / prior nodes (l21/nuclear norm dual variables)
v1_ = Composite();
weights1_ = Composite();
if init_flag
    v0_ = init_m.v0;
    weights0_ = init_m.weights0;
    s_ = Composite();
    for k = 1:K
        v1_{k} = init_m.v1(:,spectral_chunk{k});
        weights1_{k} = init_m.weights1;
        s_{k} = size(v1_{k}, 1);
    end
    fprintf('v0, v1, weigths0, weights1 uploaded \n\n')
else
    [v0_, weights0_] = hs_initialize_dual_lowrankness_serial(xsol, reweighting_alpha, sig_bar_);
    spmd
        [v1_, weights1_, s_] = hs_initialize_dual_sparsity_distributed(xsol(:,:,spectral_chunk{labindex}), Psit_, 'zpd', nlevel, reweighting_alphap, sig_);
    end
    fprintf('v0, v1, weigths0, weights1 initialized \n\n')
end

%% Data node parameters
elipse_proj_max_iter = parallel.pool.Constant(param.elipse_proj_max_iter);
elipse_proj_min_iter = parallel.pool.Constant(param.elipse_proj_min_iter);
elipse_proj_eps = parallel.pool.Constant(param.elipse_proj_eps);
adapt_eps_tol_in = parallel.pool.Constant(param.adapt_eps_tol_in);
adapt_eps_tol_out = parallel.pool.Constant(param.adapt_eps_tol_out);
adapt_eps_steps = parallel.pool.Constant(param.adapt_eps_steps);
% adapt_eps_rel_var = parallel.pool.Constant(param.adapt_eps_rel_var);
adapt_eps_change_percentage = parallel.pool.Constant(param.adapt_eps_change_percentage);

% TODO: load x directly on data nodes for Fourier transform
% ! need to be applied differently (A only exists in spmd blocks)
xi = Composite();
spmd
    Fxi_old = apply_scaled_fourier_transform(xsol(:, :, spectral_chunk{labindex}), A, No);
end

% ! TO BE UPDATED PROPERLY
if init_flag
    norm_res = Composite();
    v2_ = Composite();
    t_block = Composite();
    proj_ = Composite();
    for k = 1:K
        norm_res(k) = init_m.norm_res(k,1);
        v2_(k) = init_m.v2(k,1);
        proj_(k) = init_m.proj(k,1);
        t_block(k) = init_m.t_block(k,1);
    end
    fprintf('v2, proj, t_block, norm_res uploaded \n\n')
else
    %! this assumes the primal variable has been initialized to 0
    spmd
        [v2_, norm_res, t_block, proj_] = hs_initialize_data_worker(y);
    end
    fprintf('v2, proj, t_block, norm_res initialized \n\n')
end

% Step size for the dual variables
sigma0 = 1.0/param.nu0;
sigma1 = 1.0/param.nu1;
sigma2 = 1.0./param.nu2;

% Step size primal
tau = 0.99/3;

% Update constant dual variables
sigma00 = tau*sigma0;
sigma11 = parallel.pool.Constant(tau*sigma1);
sigma22 = Composite();
for k = 1:K
    sigma22{k} = tau*sigma2(spectral_chunk{k});
end
beta0 = param.gamma0/sigma0;
beta1 = parallel.pool.Constant(param.gamma/sigma1);
max_iter = (param.reweighting_max_iter + 1)*param.pdfb_max_iter;

% Variables for the stopping criterion
flag_convergence = 0;

if isfield(param, 'init_t_start')
    t_start = param.init_t_start;
    fprintf('t_start uploaded \n\n')
else
    param.init_t_start = 1;
    t_start = 1;
    fprintf('t_start initialized \n\n')
end

if init_flag
    rel_val = init_m.rel_val;
    end_iter = init_m.end_iter;
    t_master = init_m.t_master;
    t_l21 = init_m.t_l21;
    t_nuclear = init_m.t_nuclear;
    t_data = init_m.t_data;
    fprintf('rel_val, end_iter, t_master, t_l21, t_nuclear and t_data uploaded \n\n')
else
    rel_val = zeros(max_iter, 1);
    end_iter = zeros(max_iter, 1);
    t_master = zeros(max_iter, 1);
    t_l21 = zeros(max_iter, 1);
    t_nuclear = zeros(max_iter, 1);
    t_data = zeros(max_iter, 1);
    fprintf('rel_val, end_iter, t_master, t_l21, t_nuclear and t_data initialized \n\n')
end

if init_flag
    spmd
        [norm_residual_check_i, norm_epsilon_check_i] = sanity_check(epsilon, norm_res);
        l21_ = compute_sara_prior_distributed(xsol(:,:,spectral_chunk{labindex}), Psit_, s_);
    end  
    norm_epsilon_check = 0;
    norm_residual_check = 0;
    for k = 1:K
        norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{k};
        norm_residual_check = norm_residual_check + norm_residual_check_i{k};
    end
    norm_epsilon_check = sqrt(norm_epsilon_check);
    norm_residual_check = sqrt(norm_residual_check);

    % nuclear norm
    nuclear = nuclear_norm(xsol);

    % l21 norm
    l21 = l21_{1};

    % Log
    if (param.verbose >= 1)
        fprintf('Iter %i\n',max(t_start-1, 1));
        fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(max(t_start-1, 1)));
        fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);

        if flag_synth_data
            sol = reshape(xsol(:),numel(xsol(:))/nChannels,nChannels);
            SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
            psnrh = zeros(nChannels,1);
            for i = 1:nChannels
                psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
            end
            SNR_average = mean(psnrh);
            fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
        end
    end
end

%%
start_loop = tic;

fprintf('START THE LOOP MNRAS ver \n\n')

for t = t_start : max_iter
    
    start_iter = tic;
    
    % update primal variable
    %! to be done in parallel as well (see that afterwards)
    tw = tic;
    prev_xsol = xsol;
    xsol = max(real(xsol - g), 0);
    xhat = 2*xsol - prev_xsol;
    t_master(t) = toc(tw);

    for k = 1:K
       xi{k} = xsol(:, :, spectral_chunk{k});
    end
    
    % nuclear prior node
    tw = tic;
    [v0_, g] = hs_update_dual_lowrankness_serial(v0_, xhat, weights0_, beta0, sigma00);
    t_nuclear(t) = toc(tw);

    % update dual variables
    spmd
        % l21 prior node (full SARA prior)
        tw = tic;
        [v1_, g1_] = ...
            hs_update_dual_sparsity_distributed(v1_, Psit_, Psi_, ...
            xhat(:,:,spectral_chunk{labindex}), weights1_, ...
            beta1.Value, sigma11.Value);
        t_l21_ = toc(tw); 
        % data fidelity
        tw = tic;
        % [v2_, g2_, Fxi_old, proj_, norm_res, norm_residual_check_i, ...
        %     norm_epsilon_check_i] = ...
        %         update_dual_data_fidelity(v2_, y, xi, Fxi_old, proj_, A, At, ...
        %         G, W, pU, epsilon, elipse_proj_max_iter.Value, ...
        %         elipse_proj_min_iter.Value, elipse_proj_eps.Value, ...
        %         sigma22);
        [v2_, g2_, Fxi_old, proj_, norm_res, norm_residual_check_i, ...
            norm_epsilon_check_i] = ...
                update_dual_data_fidelity(v2_, y, xi, Fxi_old, proj_, A, At, ...
                G, W, pU, epsilon, elipse_proj_max_iter.Value, ...
                elipse_proj_min_iter.Value, elipse_proj_eps.Value, ...
                sigma22, flagDR, Sigma);
        g_ = g1_ + g2_;
        t_data_ = toc(tw); 
    end
    
    for k = 1:K
        g(:,:,spectral_chunk{k}) = g(:,:,spectral_chunk{k}) + g_{k};
    end
    
    %% Relative change of objective function
    rel_x = norm(prev_xsol(:) - xsol(:));
    norm_x = norm(xsol(:));
    rel_val(t) = rel_x/norm_x;
    end_iter(t) = toc(start_iter);
    
    % retrieve update time (average for data processes)
    t_l21(t) = 0;
    t_data(t) = 0; % just in case
    for k = 1:K
        t_l21(t) = t_l21(t) + t_l21_{k};
        t_data(t) = t_data(t) + t_data_{k};
    end
    t_data(t) = t_data(t)/K;
    t_l21(t) = t_l21(t)/K;
    
    %% Retrieve value of the monitoring variables (residual norms)
    norm_epsilon_check = 0;
    norm_residual_check = 0;
    for k = 1:K
        norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{k};
        norm_residual_check = norm_residual_check + norm_residual_check_i{k};
    end
    norm_epsilon_check = sqrt(norm_epsilon_check);
    norm_residual_check = sqrt(norm_residual_check);

    fprintf('Iter = %i, Time = %e, t_master= %e, t_l21 = %e, t_nuclear = %e, t_data = %e, rel_val = %e, epsilon = %e, residual = %e\n',t,end_iter(t),t_master(t),t_l21(t),t_nuclear(t),t_data(t),rel_val(t),norm_epsilon_check,norm_residual_check);

    %% Display
    if ~mod(t,100)

        %% compute value of the priors
        % TODO: move to l.518 if norms needs to be computed at each iteration
        nuclear = nuclear_norm(xsol);
        spmd 
            l21_ = compute_sara_prior_distributed(xsol(:,:,spectral_chunk{labindex}), Psit_, s_);
        end
        l21 = l21_{1};
        % previous_obj = obj;
        % obj = (param.gamma*l21 + param.gamma0*nuclear);
        % rel_obj = abs(previous_obj - obj)/previous_obj;
        
        % Log
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
            fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);

            if flag_synth_data
                sol = reshape(xsol(:),numel(xsol(:))/nChannels,nChannels);
                SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
                psnrh = zeros(nChannels,1);
                for i = 1:nChannels
                    psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
                end
                SNR_average = mean(psnrh);
                fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
            end
        end 
    end
    
    %% Check convergence pdfb (inner solver)
    pdfb_converged = (t - reweight_last_step_iter >= param.pdfb_min_iter) && ... % minimum number of pdfb iterations 
        ( t - reweight_last_step_iter >= param.pdfb_max_iter || ... % maximum number of pdfb iterations reached
            (rel_val(t) <= param.pdfb_rel_var && norm_residual_check <= param.pdfb_fidelity_tolerance*norm_epsilon_check) || ... % relative variation solution, objective and data fidelity within tolerance
            rel_val(t) <= param.pdfb_rel_var_low ... % relative variation really small and data fidelity criterion not satisfied yet
        );

    %% Update epsilons (in parallel)
    flag_epsilonUpdate = param.use_adapt_eps && ...  % activate espilon update 
    (t > param.adapt_eps_start) && ...               % update allowed after a minimum of iterations in the 1st reweighting
    (rel_val(t) < param.adapt_eps_rel_var);          % relative variation between 2 consecutive pdfb iterations

    if flag_epsilonUpdate
        spmd
            [epsilon, t_block] = update_epsilon(epsilon, t, t_block, ...
                norm_res, adapt_eps_tol_in.Value, ...
                adapt_eps_tol_out.Value, adapt_eps_steps.Value, ...
                adapt_eps_change_percentage.Value);
        end
    end
    %! --
    
    %% Reweighting (in parallel)
    if pdfb_converged 
        rel_x_reweighting = norm(xlast_reweight(:) - xsol(:))/norm(xlast_reweight(:));
        xlast_reweight = xsol;        
        
        reweighting_converged = pdfb_converged && ...                  % do not exit solver before the current pdfb algorithm converged
            reweight_step_count >= param.reweighting_min_iter && ...   % minimum number of reweighting iterations
            ( reweight_step_count >= param.reweighting_max_iter || ... % maximum number of reweighting iterations reached  
            rel_x_reweighting <= param.reweighting_rel_var ...         % relative variation
            );

        if reweighting_converged
            flag_convergence = 1;
            break;
        end

        fprintf('Reweighting: %i, relative variation: %e, reweighting parameter: %e \n\n', reweight_step_count+1, rel_x_reweighting, reweighting_alpha);

        % update weights nuclear norm
        weights0_ = hs_update_weights_lowrankness_serial(xsol, reweighting_alpha, sig_bar_);

        spmd
            % update weights l21-norm
            weights1_ = hs_update_weights_sparsity_distributed(xsol(:,:,spectral_chunk{labindex}), Psit_, weights1_, reweighting_alpha, sig_);

            % compute residual image
            res_ = compute_residual_images(xsol(:,:,spectral_chunk{labindex}), y, A, At, G, W, flagDR, Sigma);
        end

        %! -- TO BE CHECKED
        if flag_homotopy
            reweighting_alpha = max(param.reweighting_alpha_ff * reweighting_alpha, 1);
        end
        %! --
        param.reweighting_alpha = reweighting_alpha;
        param.init_reweight_step_count = reweight_step_count+1;
        param.init_reweight_last_iter_step = t;
        param.init_t_start = t+1;

        nuclear = nuclear_norm(xsol);
        spmd 
            l21_ = compute_sara_prior_distributed(xsol(:,:,spectral_chunk{labindex}), Psit_, s_);
        end
        l21 = l21_{1};
        % previous_obj = obj;
        % obj = (param.gamma*l21 + param.gamma0*nuclear);
        % rel_obj = abs(previous_obj - obj)/previous_obj;
        
        % compute SNR
        if flag_synth_data
            sol = reshape(xsol(:),numel(xsol(:))/nChannels,nChannels);
            SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
            psnrh = zeros(nChannels,1);
            for i = 1:nChannels
                psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
            end
            SNR_average = mean(psnrh);
            fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
        end


        if (reweight_step_count == 0) || (reweight_step_count == 1) || (~mod(reweight_step_count, param.backup_frequency))
            % Save parameters (matfile solution)
            m = matfile([name_checkpoint, '_rw=' num2str(reweight_step_count) '.mat'], ...
                'Writable', true);
            m.param = param;
            m.res = zeros(size(xsol));
            m.g = g;
            m.xsol = xsol;
            m.epsilon = cell(K, 1);
            m.v2 = cell(K, 1);
            m.proj = cell(K, 1);
            m.t_block = cell(K, 1);
            m.norm_res = cell(K, 1);
            m.v0 = v0_;
            m.v1 = zeros(s_{1}, nChannels);
            m.weights0 = weights0_;
            m.weights1 = weights1_{1};
            % Retrieve variables from workers
            % data nodes
            for k = 1:K
                m.res(:, :, spectral_chunk{k}) = res_{k};
                res_{k} = [];
                m.v2(k, 1) = v2_(k);
                m.proj(k, 1) = proj_(k);
                m.t_block(k, 1) = t_block(k);
                m.epsilon(k, 1) = epsilon(k);
                m.norm_res(k, 1) = norm_res(k);
                m.v1(:, spectral_chunk{k}) = v1_{k};
            end
            m.end_iter = end_iter;
            m.t_master = t_master;
            m.t_l21 = t_l21;
            m.t_nuclear = t_nuclear;
            m.t_data = t_data;
            m.rel_val = rel_val;
            fitswrite(m.xsol, [name_checkpoint '_xsol' '.fits'])
            fitswrite(m.res, [name_checkpoint '_res' '.fits'])
            if flag_synth_data
                m.SNR = SNR;
                m.SNR_average = SNR_average;
            end
            clear m
            
            % Log
            if (param.verbose >= 1)
                fprintf('Backup iter: %i\n',t);
                fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
                fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
                fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
            end
        end 

        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;   
        if (reweight_step_count >= param.reweighting_max_iter)
            fprintf('\n\n No more reweights \n\n');
        end     
    end
end
toc(start_loop)

% Calculate residual images
res = zeros(size(xsol));
spmd
    res_ = compute_residual_images(xsol(:,:,spectral_chunk{labindex}), y, A, At, G, W, flagDR, Sigma);
end

m = matfile([name_checkpoint, '_rw=' num2str(reweight_step_count) '.mat'], ...
    'Writable', true);
m.param = param;
m.res = zeros(size(xsol));
m.g = g;
m.epsilon = cell(K, 1);
m.v2 = cell(K, 1);
m.proj = cell(K, 1);
m.t_block = cell(K, 1);
m.norm_res = cell(K, 1);
m.v0 = v0_;
m.weights0 = weights0_;
m.v1 = zeros(s_{1}, nChannels);
m.weights1 = weights1_{1};

% Retrieve variables from workers
% data nodes
for k = 1:K
    res(:, :, spectral_chunk{k}) = res_{k};
    m.res(:, :, spectral_chunk{k}) = res_{k};
    res_{k} = [];
    m.v2(k, 1) = v2_(k);
    m.v1(:, spectral_chunk{k}) = v1_{k};
    m.proj(k, 1) = proj_(k);
    m.t_block(k, 1) = t_block(k);
    m.epsilon(k, 1) = epsilon(k);
    m.norm_res(k, 1) = norm_res(k);
end
epsilon = m.epsilon;
m.xsol = xsol;
norm_res_out = sqrt(sum(sum(sum((m.res).^2))));

% Update param structure and save
param.reweighting_alpha = reweighting_alpha;
param.init_reweight_step_count = reweight_step_count;
param.init_reweight_last_iter_step = t;
param.init_t_start = t+1;
m.param = param;
m.end_iter = end_iter;
m.t_master = t_master;
m.t_l21 = t_l21;
m.t_nuclear = t_nuclear;
m.t_data = t_data;
m.rel_val = rel_val;
fitswrite(m.xsol, [name_checkpoint '_xsol' '.fits'])
fitswrite(m.res, [name_checkpoint '_res' '.fits'])
if flag_synth_data
    sol = reshape(xsol(:),numel(xsol(:))/nChannels,nChannels);
    SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
    psnrh = zeros(nChannels,1);
    for i = 1:nChannels
        psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
    end
    SNR_average = mean(psnrh);
    m.SNR = SNR;
    m.SNR_average = SNR_average;
end
clear m

% Final log
nuclear = nuclear_norm(xsol);
spmd 
    l21_ = compute_sara_prior_distributed(xsol(:,:,spectral_chunk{labindex}), Psit_, s_);
end
l21 = l21_{1};
if (param.verbose > 0)
    if (flag_convergence == 1)
        fprintf('Solution found\n');
    else
        fprintf('Maximum number of iterations reached\n');
    end
    fprintf('Iter %i\n',t);
    fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
    fprintf('epsilon = %e, residual = %e\n', norm_epsilon_check,norm_residual_check);
    if flag_synth_data
        fprintf('SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
    end
end

end
