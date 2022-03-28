function xsol = hyperSARA(y, epsilon, A, At, pU, G, W, param, K, ...
    wavelet, nlevel, spectral_chunk, n_channels, M, N, oy, ox, ...
    warmstart_name, checkpoint_name, flag_dimensionality_reduction, Sigma, varargin)
% Implementation of the HyperSARA imaging algorithm.
%
% Implementation of the reweighting / PDFB algorithm underlying the
% HyperSARA imaging approach :cite:p:`Abdulaziz2019`. At each iteration of
% the outer reweighting scheme, the PDFB solves a problem of the form
%
% .. math::
%
%   \min_{X \in \mathbb{R}_+^{N \times L}}
%   \overline{\mu} \| X \|_{\overline{\omega}, *}
%   + \mu \|\Psi^\dagger X \|_{\omega, 2,1} + \sum_{\ell, b}
%   \iota_{ \{\| y_{\ell, b} - \Phi_{\ell, b}(\cdot) \|_2
%   \leq \varepsilon_{\ell, b}\} } (x_\ell).
%
% Parameters
% ----------
% y : composite, cell
%     Blocks of visibilities {L}{nblocks_l}.
% epsilon : composite, cell
%     Data fidelity constraints {L}{nblocks_l}.
% A : function handle
%     Measurement operator.
% At : function handle
%     Adjoint of the measurement operator.
% pU : composite, cell
%     Preconditioning matrices {L}{nblocks_l}.
% G : composite, cell of sparse matrices
%     Degridding matrices {L}{nblocks_l}.
% W : composite, cell of int[:]
%     Masks for selection of the blocks of visibilities.
% param : struct
%     Algorithm parameters.
% K : int
%     Number of data workers.
% wavelet : cell of string
%     List of wavelet basis (for the SARA prior). If the Dirac basis
%     (``'self'``) is considered, it needs to appear in last position.
% nlevel : int
%     Number of wavelet decomposition levels.
% spectral_chunk : cell of int[:]
%     List of channels handled by each data process {K, 1}.
% n_channels : int
%     Total number of spectral channels.
% M : int
%     Spatial image size along axis y.
% N : int
%     Spatial image size along axis x.
% oy : int
%     Fourier oversampling factor along axis y.
% ox : int
%     Fourier oversampling factor along axis x.
% warmstart_name : string
%     Name of a valid ``.mat`` file to initialize the solver (warm-start).
% checkpoint_name : string
%     Name defining the name for checkpoint files saved
%     throughout the iterations.
% flag_dimensionality_reduction : bool
%     Flag to indicate data dimensiomality reduction is used.
% Sigma : composite, cell of complex[:, :]
%     Dimensionality reduction matrix. (?)
% varargin{1} : double[:, :]
%     Initial value for the primal variable [M*N, L].
% varargin{2} : double[:, :]
%     Ground truth wideband image (only when debugging synthetic
%     data experiments) [M*N, L].
%
% Returns
% -------
% xsol : double[:, :]
%     Reconstructed wideband image [M*N, L].
%
%
% Important
% ---------
% - The `param` structure should contain the following fields.
%
% param.verbose (string)
%     Print log or not
%
% param.nu0 (double)
%     Norm of the splitting operator used to defined the low-rankness dual
%     variable (identity, thus fixed to 1).
% param.nu1 (double)
%     Upper bound on the norm of the SARA operator :math:`\Psi^\dagger`.
% param.nu2 (double)
%     Upper bound on the norm of the measurement operator :math:`\Phi`
% param.gamma  (double)
%     Regularization parameter (sparsity prior).
% param.gamma0  (double)
%     Regularization parameter (low-rankness prior).
%
% param.pdfb_min_iter (int)
%     Minimum number of iterations (PDFB).
% param.pdfb_max_iter (int)
%     Maximum number of iterations (PDFB).
% param.pdfb_rel_var (double)
%     Relative variation tolerance (PDFB).
% param.pdfb_fidelity_tolerance (double)
%     Tolerance to check data constraints are satisfied (PDFB).
%
% param.reweighting_max_iter (int)
%     Maximum number of reweighting steps.
% param.reweighting_min_iter (int)
%     Minimum number of reweighting steps (to reach "noise" level).
% param.reweighting_rel_var (double)
%     Tolerance relative variation (reweighting).
% param.reweighting_alpha (double)
%     Starting reweighting parameter.
% param.reweighting_sig (double)
%     Noise level (in wavelet space)
% param.reweighting_sig_bar (double)
%     Noise level (singular value space)
%
% param.elipse_proj_max_iter (int)
%     Maximum number of iterations for the FB algo that implements the
%     projection onto the ellipsoid.
% param.elipse_proj_min_iter (int)
%     Minimum number of iterations for the FB algo that implements the
%     projection onto the ellipsoid.
% parma.elipse_proj_eps (double)
%     Stopping criterion for the projection.
%
% - The following fileds are added to `param` in the course of the
%   algorithm to be able to restart it from a previous state
%
% param.reweighting_alpha (double)
%     Current state of the reweighting parameter.
% param.init_reweight_step_count (int)
%     Iteration index of the current reweighting step when the checkpoint
%     file is written.
% param.init_reweight_last_iter_step (int)
%     Global iteration index (unrolled over all pdfb iterations performed
%     in the reweighting scheme) from which to restart the algorithm.
% param.init_t_start (int)
%     Global iteration index (unrolled over all pdfb iterations performed
%     in the reweighting scheme) from which to restart the algorithm
%     (``param.init_t_start = param.init_reweight_last_iter_step + 1``).
%
% - The following fileds are added to `param` in the course of the
%   algorithm to be able to restart it from a previous state
%
% param.reweighting_alpha (double)
%     Current state of the reweighting parameter.
% param.init_reweight_step_count (int)
%     Iteration index of the current reweighting step when the checkpoint
%     file is written.
% param.init_reweight_last_iter_step (int)
%     Global iteration index (unrolled over all pdfb iterations performed
%     in the reweighting scheme) from which to restart the algorithm.
% param.init_t_start (int)
%     Global iteration index (unrolled over all pdfb iterations performed
%     in the reweighting scheme) from which to restart the algorithm
%     (``param.init_t_start = param.init_reweight_last_iter_step + 1``).
%
% Note
% ----
% The checkpoint file saved throughout the iterations is composed of the
% following variables (to be verified).
%
% param (struct)
%     Current state of the parameter structure (updated with auxiliary
%     parameters describing the state of the solver when the checkpoint
%     file has been written).
% xsol (double[:, :])
%     Current state of the image estimate.
% res (double[:, :, :])
%     Current state of the wideband residual image.
% g (double[:, :])
%     Auxiliary variable involved in the update of the primal variable.
% epsilon (cell of cell of double)
%     (Updated) Value of the l2-ball radii.
% t_block (cell of cell of int)
%     Index of the last iteration where the weigths have been updated.
% proj (cell of cell of complex[:])
%     Auxiliary variable involved in the update of the data fidelity dual
%     variables.
% norm_res (cell of cell of double)
%     Norm of the residual image (per channel per block).
% v0 (cell of double[:])
%     Dual variables associated with the low-rankness.
% v1 (double[:, :])
%     Dual variables associated with the sparsity prior.
% v2 (cell of cell complex[:])
%     Dual variables associated with the data fidelity terms (each cell
%     corresponding to a data block).
% weights0 (cell of double[:])
%     Weights associated with the low-rankness prior.
% weights1 (cell of double[:])
%     Weights associated with the sparsity prior.
% end_iter (int)
%     Last iteration (unrolled over all pdfb iterations).
% t_l21 (double)
%     Time to update the sparsity dual variable.
% t_nuclear (double)
%     Time to update the low-rankness dual variable.
% t_master (double)
%     Time to perform update/computations on the master process.
% t_data (double)
%     Time to update the data fidelity dual variables.
% rel_val (double)
%     Relative variation of the solution across the iterations.
%

%
% Deprecated fields (adaptive espilon scheme)
%
% param.reweighting_alpha_ff (double)
%     Multiplicative parameter update (< 1).
% param.use_adapt_eps (bool)
%     Flag to activate adaptive epsilon (note that there is no need to use
%     the adaptive strategy for experiments on synthetic data).
% param.adapt_eps_start (int)
%     Minimum number of iterations before starting to adjust the data
%     constaints.
% param.adapt_eps_tol_in (double)
%     Tolerance inside the :math:`\ell_2` ball (< 1).
% param.adapt_eps_tol_out (double)
%     Tolerance inside the :math:`\ell_2` ball (> 1).
% param.adapt_eps_steps (int)
%     Minimum number of iterations between consecutive data constraint
%     updates.
% param.adapt_eps_rel_var (double)
%     Bound on the relative change of the solution to trigger the update of
%     the data constraints :cite:p:`Dabbech2018`.
% param.adapt_eps_change_percentage (double)
%     Update parameter to update data constaints :cite:p:`Dabbech2018`.
%

% ------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin, A. Abdulaziz, M. Jiang
% ------------------------------------------------------------------------%
%%
% Note:
% Code based on the HyperSARA code developed by A. Abdulaziz, available at
% https://basp-group.github.io/Hyper-SARA/
% ------------------------------------------------------------------------%
%%
% ! SPMD version: use spmd for all the priors, deal with the data fidelity
% term in a single place.
%
%% REMARKS
% 1. Do we need to have the number of nodes Q, I available as parallel
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
No = M * oy * N * ox;

% instantiate Psi, Psit
spmd
    [Psi_, Psit_] = op_sp_wlt_basis(wavelet, nlevel, M, N);
end

% Initializations
init_flag = isfile(warmstart_name);
if init_flag
    init_m = matfile(warmstart_name);
    fprintf('Resume from file %s\n\n', warmstart_name);
end

% ! -- TO BE CHECKED (primal initialization)
if init_flag
    xsol = init_m.xsol;
    pdfb_rel_var_low = param.pdfb_rel_var_low;
    param = init_m.param;

    if ~isfield(param, 'pdfb_rel_var_low')
        param.pdfb_rel_var_low = pdfb_rel_var_low;
    end

    % ! -- HOMOTOPY (deactivated)
    % check flag homotopy strategy
    % if ~isfield(param, 'flag_homotopy')
    %     flag_homotopy = false;
    % else
    %     flag_homotopy = param.flag_homotopy;
    % end
    % ! --

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
    fprintf('xsol, param and epsilon uploaded \n\n');

else
    if ~isempty(varargin)
        if ~isempty(varargin{1})
            xsol = varargin{1};
        else
            xsol = zeros(M, N, n_channels);
        end

        if numel(varargin) > 1
            flag_synth_data = true;
            X0 = varargin{2};
        else
            flag_synth_data = false;
        end
    else
        xsol = zeros(M, N, n_channels);
    end
    fprintf('xsol initialized \n\n');
end
% ! assumes backup file exactly saved at the time a reweighting step occurred
xlast_reweight = xsol;
% ! --

if init_flag
    g = init_m.g;
    fprintf('g uploaded \n\n');
else
    g = zeros(size(xsol));
    fprintf('g initialized \n\n');
end

% ! Reweighting parameters
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

if isfield(param, 'init_reweight_step_count')
    reweight_step_count = param.init_reweight_step_count;
    fprintf('reweight_step_count uploaded\n\n');
else
    reweight_step_count = 0;
    fprintf('reweight_step_count initialized \n\n');
end

if isfield(param, 'init_reweight_last_iter_step')
    reweight_last_step_iter = param.init_reweight_last_iter_step;
    fprintf('reweight_last_iter_step uploaded \n\n');
else
    param.init_reweight_last_iter_step = 0;
    reweight_last_step_iter = 0;
    fprintf('reweight_last_iter_step initialized \n\n');
end
% ! --

% Primal / prior nodes (l21/nuclear norm dual variables)
v1_ = Composite();
weights1_ = Composite();
if init_flag
    v0_ = init_m.v0;
    weights0_ = init_m.weights0;
    s_ = Composite();
    for k = 1:K
        v1_{k} = init_m.v1(:, spectral_chunk{k});
        weights1_{k} = init_m.weights1;
        s_{k} = size(v1_{k}, 1);
    end
    fprintf('v0, v1, weigths0, weights1 uploaded \n\n');
else
    [v0_, weights0_] = hs_initialize_dual_lowrankness_serial(xsol, reweighting_alpha, sig_bar_);
    spmd
        [v1_, weights1_, s_] = hs_initialize_dual_sparsity_distributed(xsol(:, :, spectral_chunk{labindex}), Psit_, 'zpd', nlevel, reweighting_alphap, sig_);
    end
    fprintf('v0, v1, weigths0, weights1 initialized \n\n');
end

%% Data node parameters
elipse_proj_max_iter = parallel.pool.Constant(param.elipse_proj_max_iter);
elipse_proj_min_iter = parallel.pool.Constant(param.elipse_proj_min_iter);
elipse_proj_eps = parallel.pool.Constant(param.elipse_proj_eps);
% adapt_eps_tol_in = parallel.pool.Constant(param.adapt_eps_tol_in);
% adapt_eps_tol_out = parallel.pool.Constant(param.adapt_eps_tol_out);
% adapt_eps_steps = parallel.pool.Constant(param.adapt_eps_steps);
% adapt_eps_rel_var = parallel.pool.Constant(param.adapt_eps_rel_var);
% adapt_eps_change_percentage = parallel.pool.Constant(param.adapt_eps_change_percentage);

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
        norm_res(k) = init_m.norm_res(k, 1);
        v2_(k) = init_m.v2(k, 1);
        proj_(k) = init_m.proj(k, 1);
        t_block(k) = init_m.t_block(k, 1);
    end
    fprintf('v2, proj, t_block, norm_res uploaded \n\n');
else
    % ! this assumes the primal variable has been initialized to 0
    spmd
        [v2_, norm_res, t_block, proj_] = hs_initialize_data_worker(y);
    end
    fprintf('v2, proj, t_block, norm_res initialized \n\n');
end

% Step size for the dual variables
sigma0 = 1.0 / param.nu0;
sigma1 = 1.0 / param.nu1;
sigma2 = 1.0 ./ param.nu2;

% Step size primal
tau = 0.99 / 3;

% Update constant dual variables
sigma00 = tau * sigma0;
sigma11 = parallel.pool.Constant(tau * sigma1);
sigma22 = Composite();
for k = 1:K
    sigma22{k} = tau * sigma2(spectral_chunk{k});
end
beta0 = param.gamma0 / sigma0;
beta1 = parallel.pool.Constant(param.gamma / sigma1);
max_iter = (param.reweighting_max_iter + 1) * param.pdfb_max_iter;

% Variables for the stopping criterion
flag_convergence = 0;

if isfield(param, 'init_t_start')
    t_start = param.init_t_start;
    fprintf('t_start uploaded \n\n');
else
    param.init_t_start = 1;
    t_start = 1;
    fprintf('t_start initialized \n\n');
end

if init_flag
    rel_val = init_m.rel_val;
    end_iter = init_m.end_iter;
    t_master = init_m.t_master;
    t_l21 = init_m.t_l21;
    t_nuclear = init_m.t_nuclear;
    t_data = init_m.t_data;
    fprintf('rel_val, end_iter, t_master, t_l21, t_nuclear and t_data uploaded \n\n');
else
    rel_val = zeros(max_iter, 1);
    end_iter = zeros(max_iter, 1);
    t_master = zeros(max_iter, 1);
    t_l21 = zeros(max_iter, 1);
    t_nuclear = zeros(max_iter, 1);
    t_data = zeros(max_iter, 1);
    fprintf('rel_val, end_iter, t_master, t_l21, t_nuclear and t_data initialized \n\n');
end

if init_flag
    spmd
        [norm_residual_check_i, norm_epsilon_check_i] = sanity_check(epsilon, norm_res);
        l21_ = compute_sara_prior_distributed(xsol(:, :, spectral_chunk{labindex}), Psit_, s_);
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
    if param.verbose >= 1
        fprintf('Iter %i\n', max(t_start - 1, 1));
        fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(max(t_start - 1, 1)));
        fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);

        if flag_synth_data
            sol = reshape(xsol(:), numel(xsol(:)) / n_channels, n_channels);
            SNR = 20 * log10(norm(X0(:)) / norm(X0(:) - sol(:)));
            psnrh = zeros(n_channels, 1);
            for i = 1:n_channels
                psnrh(i) = 20 * log10(norm(X0(:, i)) / norm(X0(:, i) - sol(:, i)));
            end
            SNR_average = mean(psnrh);
            fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
        end
    end
end

%%
start_loop = tic;

fprintf('START THE LOOP MNRAS ver \n\n');

for t = t_start:max_iter

    start_iter = tic;

    % update primal variable
    % ! to be done in parallel as well (see that afterwards)
    tw = tic;
    prev_xsol = xsol;
    xsol = max(real(xsol - g), 0);
    xhat = 2 * xsol - prev_xsol;
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
            xhat(:, :, spectral_chunk{labindex}), weights1_, ...
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
                sigma22, flag_dimensionality_reduction, Sigma);
        g_ = g1_ + g2_;
        t_data_ = toc(tw);
    end

    for k = 1:K
        g(:, :, spectral_chunk{k}) = g(:, :, spectral_chunk{k}) + g_{k};
    end

    %% Relative change of objective function
    rel_x = norm(prev_xsol(:) - xsol(:));
    norm_x = norm(xsol(:));
    rel_val(t) = rel_x / norm_x;
    end_iter(t) = toc(start_iter);

    % retrieve update time (average for data processes)
    t_l21(t) = 0;
    t_data(t) = 0; % just in case
    for k = 1:K
        t_l21(t) = t_l21(t) + t_l21_{k};
        t_data(t) = t_data(t) + t_data_{k};
    end
    t_data(t) = t_data(t) / K;
    t_l21(t) = t_l21(t) / K;

    %% Retrieve value of the monitoring variables (residual norms)
    norm_epsilon_check = 0;
    norm_residual_check = 0;
    for k = 1:K
        norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{k};
        norm_residual_check = norm_residual_check + norm_residual_check_i{k};
    end
    norm_epsilon_check = sqrt(norm_epsilon_check);
    norm_residual_check = sqrt(norm_residual_check);

    fprintf('Iter = %i, Time = %e, t_master= %e, t_l21 = %e, t_nuclear = %e, t_data = %e, rel_val = %e, epsilon = %e, residual = %e\n', t, end_iter(t), t_master(t), t_l21(t), t_nuclear(t), t_data(t), rel_val(t), norm_epsilon_check, norm_residual_check);

    %% Display
    if ~mod(t, 100)

        %% compute value of the priors
        % TODO: move to l.518 if norms needs to be computed at each iteration
        nuclear = nuclear_norm(xsol);
        spmd
            l21_ = compute_sara_prior_distributed(xsol(:, :, spectral_chunk{labindex}), Psit_, s_);
        end
        l21 = l21_{1};
        % previous_obj = obj;
        % obj = (param.gamma*l21 + param.gamma0*nuclear);
        % rel_obj = abs(previous_obj - obj)/previous_obj;

        % Log
        if param.verbose >= 1
            fprintf('Iter %i\n', t);
            fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
            fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);

            if flag_synth_data
                sol = reshape(xsol(:), numel(xsol(:)) / n_channels, n_channels);
                SNR = 20 * log10(norm(X0(:)) / norm(X0(:) - sol(:)));
                psnrh = zeros(n_channels, 1);
                for i = 1:n_channels
                    psnrh(i) = 20 * log10(norm(X0(:, i)) / norm(X0(:, i) - sol(:, i)));
                end
                SNR_average = mean(psnrh);
                fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
            end
        end
    end

    %% Check convergence pdfb (inner solver)
    pdfb_converged = (t - reweight_last_step_iter >= param.pdfb_min_iter) && ... % minimum number of pdfb iterations
        (t - reweight_last_step_iter >= param.pdfb_max_iter || ... % maximum number of pdfb iterations reached
            (rel_val(t) <= param.pdfb_rel_var && norm_residual_check <= param.pdfb_fidelity_tolerance * norm_epsilon_check) || ... % relative variation solution, objective and data fidelity within tolerance
            rel_val(t) <= param.pdfb_rel_var_low ... % relative variation really small and data fidelity criterion not satisfied yet
        );

    % %% Update epsilons (in parallel)
    % flag_epsilonUpdate = param.use_adapt_eps && ...  % activate espilon update
    % (t > param.adapt_eps_start) && ...               % update allowed after a minimum of iterations in the 1st reweighting
    % (rel_val(t) < param.adapt_eps_rel_var);          % relative variation between 2 consecutive pdfb iterations
    %
    % if flag_epsilonUpdate
    %     spmd
    %         [epsilon, t_block] = update_epsilon(epsilon, t, t_block, ...
    %             norm_res, adapt_eps_tol_in.Value, ...
    %             adapt_eps_tol_out.Value, adapt_eps_steps.Value, ...
    %             adapt_eps_change_percentage.Value);
    %     end
    % end
    % ! --

    %% Reweighting (in parallel)
    if pdfb_converged
        rel_x_reweighting = norm(xlast_reweight(:) - xsol(:)) / norm(xlast_reweight(:));
        xlast_reweight = xsol;

        reweighting_converged = pdfb_converged && ...                  % do not exit solver before the current pdfb algorithm converged
            reweight_step_count >= param.reweighting_min_iter && ...   % minimum number of reweighting iterations
            (reweight_step_count >= param.reweighting_max_iter || ... % maximum number of reweighting iterations reached
            rel_x_reweighting <= param.reweighting_rel_var ...         % relative variation
            );

        if reweighting_converged
            flag_convergence = 1;
            break
        end

        fprintf('Reweighting: %i, relative variation: %e, reweighting parameter: %e \n\n', reweight_step_count + 1, rel_x_reweighting, reweighting_alpha);

        % update weights nuclear norm
        weights0_ = hs_update_weights_lowrankness_serial(xsol, reweighting_alpha, sig_bar_);

        spmd
            % update weights l21-norm
            weights1_ = hs_update_weights_sparsity_distributed(xsol(:, :, spectral_chunk{labindex}), Psit_, weights1_, reweighting_alpha, sig_);

            % compute residual image
            res_ = compute_residual_images(xsol(:, :, spectral_chunk{labindex}), y, A, At, G, W, flag_dimensionality_reduction, Sigma);
        end

        % ! -- HOMOTOPY (deactivated)
        % if flag_homotopy
        %     reweighting_alpha = max(param.reweighting_alpha_ff * reweighting_alpha, 1);
        % end
        % ! --
        param.reweighting_alpha = reweighting_alpha;
        param.init_reweight_step_count = reweight_step_count + 1;
        param.init_reweight_last_iter_step = t;
        param.init_t_start = t + 1;

        nuclear = nuclear_norm(xsol);
        spmd
            l21_ = compute_sara_prior_distributed(xsol(:, :, spectral_chunk{labindex}), Psit_, s_);
        end
        l21 = l21_{1};
        % previous_obj = obj;
        % obj = (param.gamma*l21 + param.gamma0*nuclear);
        % rel_obj = abs(previous_obj - obj)/previous_obj;

        % compute SNR
        if flag_synth_data
            sol = reshape(xsol(:), numel(xsol(:)) / n_channels, n_channels);
            SNR = 20 * log10(norm(X0(:)) / norm(X0(:) - sol(:)));
            psnrh = zeros(n_channels, 1);
            for i = 1:n_channels
                psnrh(i) = 20 * log10(norm(X0(:, i)) / norm(X0(:, i) - sol(:, i)));
            end
            SNR_average = mean(psnrh);
            fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
        end

        if (reweight_step_count == 0) || (reweight_step_count == 1) || (~mod(reweight_step_count, param.backup_frequency))
            % Save parameters (matfile solution)
            m = matfile([checkpoint_name, '_rw=' num2str(reweight_step_count) '.mat'], ...
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
            m.v1 = zeros(s_{1}, n_channels);
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
            fitswrite(m.xsol, [checkpoint_name '_xsol' '.fits']);
            fitswrite(m.res, [checkpoint_name '_res' '.fits']);
            if flag_synth_data
                m.SNR = SNR;
                m.SNR_average = SNR_average;
            end
            clear m;

            % Log
            if param.verbose >= 1
                fprintf('Backup iter: %i\n', t);
                fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
                fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
                if flag_synth_data
                    fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
                end
            end
        end

        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;
        if reweight_step_count >= param.reweighting_max_iter
            fprintf('\n\n No more reweights \n\n');
        end
    end
end
toc(start_loop);

% Calculate residual images
res = zeros(size(xsol));
spmd
    res_ = compute_residual_images(xsol(:, :, spectral_chunk{labindex}), y, A, At, G, W, flag_dimensionality_reduction, Sigma);
end

m = matfile([checkpoint_name, '_rw=' num2str(reweight_step_count) '.mat'], ...
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
m.v1 = zeros(s_{1}, n_channels);
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
param.init_t_start = t + 1;
m.param = param;
m.end_iter = end_iter;
m.t_master = t_master;
m.t_l21 = t_l21;
m.t_nuclear = t_nuclear;
m.t_data = t_data;
m.rel_val = rel_val;
fitswrite(m.xsol, [checkpoint_name '_xsol' '.fits']);
fitswrite(m.res, [checkpoint_name '_res' '.fits']);
if flag_synth_data
    sol = reshape(xsol(:), numel(xsol(:)) / n_channels, n_channels);
    SNR = 20 * log10(norm(X0(:)) / norm(X0(:) - sol(:)));
    psnrh = zeros(n_channels, 1);
    for i = 1:n_channels
        psnrh(i) = 20 * log10(norm(X0(:, i)) / norm(X0(:, i) - sol(:, i)));
    end
    SNR_average = mean(psnrh);
    m.SNR = SNR;
    m.SNR_average = SNR_average;
end
clear m;

% Final log
nuclear = nuclear_norm(xsol);
spmd
    l21_ = compute_sara_prior_distributed(xsol(:, :, spectral_chunk{labindex}), Psit_, s_);
end
l21 = l21_{1};
if param.verbose > 0
    if flag_convergence == 1
        fprintf('Solution found\n');
    else
        fprintf('Maximum number of iterations reached\n');
    end
    fprintf('Iter %i\n', t);
    fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
    fprintf('epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
    if flag_synth_data
        fprintf('SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
    end
end

end
