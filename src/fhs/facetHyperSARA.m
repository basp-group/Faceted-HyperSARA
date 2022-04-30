function xsol = facetHyperSARA(y, epsilon, A, At, pU, G, W, param, ...
    Qx, Qy, K, wavelet, filter_length, nlevel, window_type, ...
    spectral_chunk, n_channels, overlap_size, ...
    M, N, oy, ox, warmstart_name, checkpoint_name, ...
    flag_dimensionality_reduction, Lambda, ...
    varargin)
% Implementation of the Faceted HyperSARA imaging algorithm.
%
% Implementation of the reweighting / PDFB algorithm underlying the
% Faceted HyperSARA imaging approach :cite:p:`Thouvenin2021`. At
% each iteration of the outer reweighting scheme, the PDFB solves a problem
% of the form
%
% .. math::
%
%   \min_{X \in \mathbb{R}_+^{N \times L}}
%   \sum_{q=1}^Q \big(
%   \overline{\mu}_q \| D_q \overline{S}_q X \|_{\overline{\omega}_q, *}
%   + \mu \|\Psi_q^\dagger S_q X \|_{\omega_q, 2,1} \big)
%   + \sum_{\ell, b} \iota_{ \{\| y_{\ell, b} - \Phi_{\ell, b}(\cdot) \|_2
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
% Qx : int
%     Number of spatial facets along spatial axis x.
% Qy : int
%     Number of spatial facets along spatial axis y.
% K : int
%     Number of data workers.
% wavelet : cell of string
%     List of wavelet basis (for the SARA prior). If the Dirac basis
%     (``'self'``) is considered, it needs to appear in last position.
% filter_length : int[:]
%     Length of each wavelet filter considered for the SARA dictionary (0
%     by convention for the Dirac basis).
% nlevel : int
%     Number of wavelet decomposition levels.
% window_type : string (``"triangular"`` by default, ``"hamming"`` or
%     ``"pc"`` (piecewise-constant))
%     Type of apodization window considered for the faceted nuclear norm
%     prior.
% spectral_chunk : cell of int[:]
%     List of indices of spectral channels handled by each of the `K` data
%     workers
% n_channels : int
%     Total number of spectral channels.
% overlap_size : int[2]
%     Overlap size between consecutive facets along each axis (y and x) for
%     the faceted low-rankness prior.
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
% Lambda : composite, cell of complex[:, :]
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
% Important
% ---------
% - The `param` structure should contain the following fields.
%
% param.verbose (string)
%     Print log or not
% param.save_intermediate_results_mat (bool)
%     Flag to save all estimates in a ``.mat`` file after each re-weighted
%     problem. The file can be used as a warmstart.
% param.nu0 (double)
%     Norm of the splitting operator used to defined the faceted
%     low-rankness dual variable (identity for ``rectangular`` window, thus
%     fixed to 1).
% param.nu1 (double)
%     Upper bound on the norm of the SARA operator :math:`\Psi^\dagger`.
% param.nu2 (double)
%     Upper bound on the norm of the measurement operator :math:`\Phi`
% param.gamma  (double)
%     Regularization parameter (sparsity prior).
% param.gamma0  (double[:)
%     Regularization parameter (faceted low-rankness prior).
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
% - The following fileds are involved in the adaptive estimation of the
%   noise level
%
% param.adjust_flag_noise (bool)
%     Flag to activate the adaptive procedure to estimate the noise level.
% param.adjust_noise_min_iter (int)
%     Minimum number of iterations  to enable the adjustement of the
%     estimate of the noise level.
% param.adjust_noise_rel_var (double)
%     Tolerance relative variation (reweighting) to  to enable the adjustement of the
%     estimate of the noise level.
% param.adjust_noise_start_iter (int)
%     Number of iterations to force triggering noise adjustement.
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
%     Current state of the image estimate [M*N, L].
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
%     Dual variables associated with the low-rankness prior {Q}.
% v1 (cell of double[:, :])
%     Dual variables associated with the sparsity prior {Q}.
% v2 (cell of cell complex[:])
%     Dual variables associated with the data fidelity terms (each cell
%     corresponding to a data block) {K}.
% weights0 (cell of double[:])
%     Weights associated with the faceted low-rankness prior {Q}.
% weights1 (cell of double[:])
%     Weights associated with the sparsity prior {Q}.
% end_iter (int)
%     Last iteration (unrolled over all pdfb iterations).
% t_facet (double)
%     Time to update all the variables assoacited with a spatial facet..
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
% Version with a fixed overlap for the faceted nuclear norm, larger or
% smaller than the extension needed for the 2D segmented discrete wavelet
% transforms (sdwt2). Includes spatial weihting correction for the faceted
% nuclear norm (triangular, hamming, piecewise_constant, no weights by
% default).
% ------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin, A. Abdulaziz, A. Dabbech, M. Jiang
% ------------------------------------------------------------------------%
%%
% Note:
% Code based on the HyperSARA code developed by A. Abdulaziz, available at
% https://basp-group.github.io/Hyper-SARA/
% ------------------------------------------------------------------------%
%%
% SPMD version: use spmd for all the priors, deal with the data fidelity
% term in a  le place. Constant overlap for the nuclear norm assuming
% overlap_size smaller than the smallest overlap for the sdwt2 (the other option
% would also change the com cation process (borders and reduction
% operation)). overlap_size power(2, nlevel)-1)*(max(filter_length(:)-1))
%
%% NOTE
% this version relies on a specialised version of sdwt2, slightly less
% general but faster (based on Arwa's work).
%%

% size of the oversampled Fourier space (vectorized)
No = M * oy * N * ox;

% -- instantiate auxiliary variables for sdwt2
% define reference 2D facets (no overlap)
Q = Qx * Qy;
rg_y = split_range(Qy, M);
rg_x = split_range(Qx, N);
I = zeros(Q, 2);
dims = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx - 1) * Qy + qy;
        I(q, :) = [rg_y(qy, 1) - 1, rg_x(qx, 1) - 1];
        dims(q, :) = [rg_y(qy, 2) - rg_y(qy, 1) + 1, rg_x(qx, 2) - rg_x(qx, 1) + 1];
    end
end
clear rg_y rg_x;

rg_yo = split_range(Qy, M, overlap_size(1));
rg_xo = split_range(Qx, N, overlap_size(2));
Io = zeros(Q, 2);
dims_o = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx - 1) * Qy + qy;
        Io(q, :) = [rg_yo(qy, 1) - 1, rg_xo(qx, 1) - 1];
        dims_o(q, :) = [rg_yo(qy, 2) - rg_yo(qy, 1) + 1, rg_xo(qx, 2) - rg_xo(qx, 1) + 1];
    end
end
clear rg_yo rg_xo;

% instantiate auxiliary variables for faceted wavelet transforms involved
% in SARA (sdwt2)
[~, dims_overlap_ref, I_overlap, dims_overlap, status, offset, offsetL, ...
    offsetR, Ncoefs, temLIdxs, temRIdxs] = sdwt2_setup([M, N], I, dims, nlevel, wavelet, filter_length);

% define parallel constants (known by each worker)
Qyp = parallel.pool.Constant(Qy);
Qxp = parallel.pool.Constant(Qx);
Qp = parallel.pool.Constant(Q);
Kp = parallel.pool.Constant(K);
spectral_chunkp = parallel.pool.Constant(spectral_chunk);
waveletp = parallel.pool.Constant(wavelet);
nlevelp = parallel.pool.Constant(nlevel);
offsetp = parallel.pool.Constant(offset);

% define auxiliary composite variables (local to a given worker)
[Iq, dims_q, dims_oq, dims_overlap_ref_q, I_overlap_q, ...
    dims_overlap_q, status_q, offsetLq, offsetRq, Ncoefs_q, temLIdxs_q, ...
    temRIdxs_q, overlap_g_south, overlap_g_east, overlap_g_south_east, ...
    overlap, apodization_window, crop_low_rank, crop_sparsity] = fhs_setup_priors(Qx, Qy, ...
    I, dims, dims_o, dims_overlap_ref, I_overlap, dims_overlap, status, ...
    offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs, window_type, overlap_size);

% Initializations
init_flag = isfile(warmstart_name);
if init_flag
    init_m = matfile(warmstart_name);
    fprintf('Resume from file %s\n\n', warmstart_name);
end

% ! -- HOMOTOPY (deactivated)
% if ~isfield(param, 'flag_homotopy')
%     flag_homotopy = false;
% else
%     flag_homotopy = param.flag_homotopy;
% end
% ! --

% ! -- Primal initialization
if init_flag
    xsol = init_m.xsol;
    pdfb_rel_var_low = param.pdfb_rel_var_low;
    param = init_m.param;

    if ~isfield(param, 'pdfb_rel_var_low')
        param.pdfb_rel_var_low = pdfb_rel_var_low;
    end

    epsilon = Composite();
    for k = 1:K
        epsilon(Q + k) = init_m.epsilon(k, 1);
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
% ! --

g_q = Composite();
xsol_q = Composite();
if init_flag
    for q = 1:Q
        xsol_q{q} = xsol(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :);
        g_q{q} = init_m.g(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :);
    end
    fprintf('g uploaded \n\n');
else
    for q = 1:Q
        xsol_q{q} = xsol(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :);
        g_q{q} = zeros([dims(q, :), n_channels]);
    end
    fprintf('g initialized \n\n');
end

% ! -- TO BE CHECKED
% Reweighting parameters
sig_ = Composite();
sig_bar_ = Composite();
for q = 1:Q
    sig_bar_{q} = param.reweighting_sig_bar(q);
    sig_{q} = param.reweighting_sig;
end
reweighting_alpha = param.reweighting_alpha;
reweighting_alphap = Composite();
for q = 1:Q
    reweighting_alphap{q} = reweighting_alpha;
end
% reweighting_alpha_ffp = parallel.pool.Constant(param.reweighting_alpha_ff);
reweighting_flag = param.reweighting_flag;

if isfield(param, 'init_reweight_step_count')
    reweight_step_count = param.init_reweight_step_count;
    fprintf('reweight_step_count uploaded\n\n');
else
    param.init_reweight_step_count = 0;
    reweight_step_count = 0;
    fprintf('reweight_step_count initialized \n\n');
end

if isfield(param, 'init_reweight_last_iter_step')
    redefine_min_task_last_step_iter = param.init_reweight_last_iter_step;
    fprintf('reweight_last_iter_step uploaded \n\n');
else
    param.init_reweight_last_iter_step = 0;
    redefine_min_task_last_step_iter = 0;
    fprintf('reweight_last_iter_step initialized \n\n');
end
% ! --

% Primal / prior nodes (l21/nuclear norm dual variables)
% ! assumes backup file exactly saved at the time a reweighting step occured
% ! (initialized to xsol_q)
v0_ = Composite();
weights0_ = Composite();
v1_ = Composite();
weights1_ = Composite();
xlast_reweight_q = Composite();
if init_flag
    for q = 1:Q
        v0_(q) = init_m.v0(q, 1);
        v1_(q) = init_m.v1(q, 1);
        weights0_(q) = init_m.weights0(q, 1);
        weights1_(q) = init_m.weights1(q, 1);
    end
    spmd
        if labindex <= Qp.Value
            max_dims = max(dims_overlap_ref_q, dims_oq);
            xlast_reweight_q = xsol_q;
        end
    end
    fprintf('v0, v1, weigths0, weights1 uploaded \n\n');
else
    spmd
        if labindex <= Qp.Value
            max_dims = max(dims_overlap_ref_q, dims_oq);
            xlast_reweight_q = xsol_q;
            % !-- TO BE CHECKED
            x_facet = zeros([max_dims, size(xsol_q, 3)]);
            x_facet(overlap(1) + 1:end, overlap(2) + 1:end, :) = xsol_q;
            x_facet = comm2d_update_borders(x_facet, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);

            % weights initialized from initial primal variable, dual variables to 0
            [v0_, v1_, weights0_, weights1_] = fhs_initialize_dual_and_weights(x_facet, ...
                Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, max_dims - crop_low_rank, n_channels, dims_overlap_ref_q, ...
                offsetLq, offsetRq, reweighting_alphap, crop_sparsity, crop_low_rank, apodization_window, sig_, sig_bar_);

            % weights and dual variables initialized from initial primal variable
            % [v0, v1, weights0, weights1] = fhs_initialize_dual and_weights2(x_facet, ...
            %     Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
            %     offsetLq, offsetRq, reweighting_alphap, crop_sparsity, crop_low_rank, apodization_window);
            % !--
        end
    end
    fprintf('v0, v1, weigths0, weights1 initialized \n\n');
end

%% Adptive noise level estimation
if isfield(param, 'adjust_flag_noise')
    adjust_flag_noise = param.adjust_flag_noise;
else; adjust_flag_noise = false;
end
if isfield(param, 'adjust_noise_min_iter')
    adjust_noise_min_iter = param.adjust_noise_min_iter;
else; adjust_noise_min_iter = 100;
end
if isfield(param, 'adjust_noise_rel_var')
    adjust_noise_rel_var = param.adjust_noise_rel_var;
else; adjust_noise_rel_var = 1e-3;
end
if isfield(param, 'adjust_noise_start_iter')
    adjust_noise_start_iter = param.adjust_noise_start_iter;
else; adjust_noise_start_iter = 2000;
end
if isfield(param, 'adjust_noise_change_percentage')
    adjust_noise_change_percentage = param.adjust_noise_change_percentage;
else; adjust_noise_change_percentage = 0.5;
end

if isfield(param, 'adjust_noise_start_change_percentage')
    adjust_noise_start_change_percentage = param.adjust_noise_start_change_percentage;
else; adjust_noise_start_change_percentage = 0.1;
end
%% Data node parameters
elipse_proj_max_iter = parallel.pool.Constant(param.elipse_proj_max_iter);
elipse_proj_min_iter = parallel.pool.Constant(param.elipse_proj_min_iter);
elipse_proj_eps = parallel.pool.Constant(param.elipse_proj_eps);

% adapt_eps_tol_in = parallel.pool.Constant(param.adapt_eps_tol_in);
% adapt_eps_tol_out = parallel.pool.Constant(param.adapt_eps_tol_out);
% adapt_eps_steps = parallel.pool.Constant(param.adapt_eps_steps);
% adapt_eps_change_percentage = parallel.pool.Constant(param.adapt_eps_change_percentage);

if init_flag
    norm_res = Composite();
    v2_ = Composite();
    t_block = Composite();
    proj_ = Composite();
    for k = 1:K
        norm_res(Q + k) = init_m.norm_res(k, 1);
        v2_(Q + k) = init_m.v2(k, 1);
        proj_(Q + k) = init_m.proj(k, 1);
        t_block(Q + k) = init_m.t_block(k, 1);
    end
    fprintf('v2, proj, t_block, norm_res uploaded \n\n');
else
    % ! assumes primal variable initialized to 0
    spmd
        if labindex > Qp.Value
            [v2_, norm_res, t_block, proj_] = fhs_initialize_data_worker(y);
        end
    end
    fprintf('v2, proj, t_block, norm_res initialized \n\n');
end

% initialize xi and Fxi_old
spmd
    if labindex <= Qp.Value
        % send xhat_q (communication towards the data nodes)
        for i = 1:K
            labSend(xsol_q(:, :, spectral_chunkp.Value{i}), Qp.Value + i);
        end
    else
        xi = zeros(M, N, numel(y));
        for q = 1:Qp.Value
            xi(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :) = ...
                labReceive(q);
        end

        Fxi_old = zeros(No, numel(y));
        for l = 1:numel(y)
            Fxi_old(:, l) = A(xi(:, :, l));
        end
    end
end

% Step sizes for the dual variables
sigma0 = 1.0 / param.nu0;
sigma1 = 1.0 / param.nu1;
sigma2 = 1.0 ./ param.nu2;

% Step size primal
tau = 0.99 / 3;

% Update constant dual variables
sigma00 = parallel.pool.Constant(tau * sigma0);
sigma11 = parallel.pool.Constant(tau * sigma1);
% sigma22 = parallel.pool.Constant(tau*sigma2);
sigma22 = Composite();
for k = 1:K
    sigma22{Q + k} = tau * sigma2(spectral_chunk{k});
end
% beta0 = parallel.pool.Constant(param.gamma0/sigma0);
beta0 = Composite();
for q = 1:Q
    beta0{q} = param.gamma0(q) / sigma0;
end
% % beta1 = parallel.pool.Constant(param.gamma / sigma1);
beta1 = (param.gamma / sigma1);

% already in param!
% param.alph = alph;
% param.alph_bar = alph_bar;
alph = param.alph;
alph_bar = param.alph_bar;
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

max_iter = (param.reweighting_max_iter + 1) * param.pdfb_max_iter;
if init_flag
    rel_val = init_m.rel_val;
    end_iter = init_m.end_iter;
    t_facet = init_m.t_facet;
    t_data = init_m.t_data;
    fprintf('rel_val, end_iter, t_facet and t_data uploaded \n\n');
else
    rel_val = zeros(max_iter, 1);
    end_iter = zeros(max_iter, 1);
    t_facet = zeros(max_iter, 1);
    t_data = zeros(max_iter, 1);
    fprintf('rel_val, end_iter, t_facet and t_data initialized \n\n');
end

% ! check warm-start worked as expected
if init_flag
    spmd
        if labindex > Qp.Value
            [norm_residual_check_i, norm_epsilon_check_i] = sanity_check(epsilon, norm_res);
        end
    end
    norm_epsilon_check = 0;
    norm_residual_check = 0;
    for i = Q + 1:Q + K
        norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{i};
        norm_residual_check = norm_residual_check + norm_residual_check_i{i};
    end
    norm_epsilon_check = sqrt(norm_epsilon_check);
    norm_residual_check = sqrt(norm_residual_check);

    % compute value of the priors in parallel
    spmd
        if labindex <= Qp.Value
            % compute values for the prior terms
            x_facet = zeros([max_dims, size(xsol_q, 3)]);
            x_facet(overlap(1) + 1:end, overlap(2) + 1:end, :) = xsol_q;
            x_facet = comm2d_update_borders(x_facet, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
            [l21_norm, nuclear_norm] = fhs_compute_facet_prior(x_facet, Iq, ...
                offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
                offsetLq, offsetRq, crop_sparsity, crop_low_rank, apodization_window, size(v1_));
        end
    end

    % retrieve value of the priors
    l21 = 0;
    nuclear = 0;
    for q = 1:Q
        l21 = l21 + l21_norm{q};
        nuclear = nuclear + nuclear_norm{q};
    end
    % obj = param.gamma0*nuclear + param.gamma*l21;

    % Log
    if param.verbose >= 1
        fprintf('Iter %i\n', max(t_start - 1, 1));
        fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(max(t_start - 1, 1)));
        fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);

        if flag_synth_data
            for q = 1:Q
                xsol(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :) = xsol_q{q};
            end
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

start_loop = tic;

fprintf('START THE LOOP MNRAS ver \n\n');

for t = t_start:max_iter

    start_iter = tic;
    if  adjust_flag_noise && t == 1; beta1 = 0; % de-activate prior
    end
    spmd
        if adjust_flag_noise && t == 1; beta0 = 0; % de-activate prior
        end
        if labindex <= Qp.Value
            % --- primal/prior nodes (1:Q)

            % update primal variable
            tw = tic;
            [xsol_q, xhat_q, rel_x_q, norm_x_q] = fhs_update_primal(xsol_q, g_q);
            g_q = []; % mem reasons
            t_op = toc(tw);

            % send xhat_q (communication towards the data nodes)
            for k = 1:K
                labSend(xsol_q(:, :, spectral_chunkp.Value{k}), Qp.Value + k);
            end

            % update borders (-> versions of xhat with overlap)
            tw = tic;
            x_facet = zeros([max_dims, size(xsol_q, 3)]);
            if (beta1 + beta0) > 0
                x_facet(overlap(1) + 1:end, overlap(2) + 1:end, :) = xhat_q;  xhat_q = []; % mem reasons
                x_facet = comm2d_update_borders(x_facet, overlap, overlap_g_south_east, ...
                    overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
            end
            % update dual variables (nuclear, l21)
            if beta0 > 0 % check if prior is active
                [v0_, g0] = fhs_update_dual_lowrankness(v0_, ...
                    x_facet(crop_low_rank(1) + 1:end, crop_low_rank(2) + 1:end, :), ...
                    apodization_window, weights0_, beta0);
            end
            if beta1 > 0 % check if prior is active
                [v1_, g1] = fhs_update_dual_sparsity(v1_, ...
                    x_facet(crop_sparsity(1) + 1:end, crop_sparsity(2) + 1:end, :), ...
                    weights1_, beta1, Iq, ...
                    dims_q, I_overlap_q, dims_overlap_q, offsetp.Value, status_q, ...
                    nlevelp.Value, waveletp.Value, Ncoefs_q, temLIdxs_q, ...
                    temRIdxs_q, offsetLq, offsetRq, dims_overlap_ref_q);
            end

            g = zeros(size(x_facet));
            if beta0 > 0 % check if prior is active
                g(crop_low_rank(1) + 1:end, crop_low_rank(2) + 1:end, :) = sigma00.Value * g0;   g0 = []; % mem reasons
            end
            if beta1 > 0 % check if prior is active
                g(crop_sparsity(1) + 1:end, crop_sparsity(2) + 1:end, :) = sigma11.Value * g1 + ...
                    g(crop_sparsity(1) + 1:end, crop_sparsity(2) + 1:end, :);  g1 = []; % mem reasons
            end
            g = comm2d_reduce(g, overlap, Qyp.Value, Qxp.Value);
            t_op = t_op + toc(tw);

            % compute g_ for the final update term
            g_q = g(overlap(1) + 1:end, overlap(2) + 1:end, :);  g = []; % mem reasons

            % retrieve portions of g from the data nodes
            for k = 1:Kp.Value
                g_q(:, :, spectral_chunkp.Value{k}) = g_q(:, :, spectral_chunkp.Value{k}) + ...
                    labReceive(Qp.Value + k);
            end
        else
            % data nodes (Q+1:Q+K) (no parallelisation over data blocks, just frequency)
            % retrieve xhat_i from the prior/primal nodes
            for q = 1:Qp.Value
                xi(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :) = labReceive(q);
            end
            tw = tic;
            [v2_, g2, Fxi_old, proj_, norm_res, norm_residual_check_i, norm_epsilon_check_i] = update_dual_data_fidelity(v2_, y, xi, ...
                Fxi_old, proj_, A, At, G, W, pU, epsilon, ...
                elipse_proj_max_iter.Value, elipse_proj_min_iter.Value, ...
                elipse_proj_eps.Value, sigma22, flag_dimensionality_reduction, Lambda); %#ok<*ASGLU>
            t_op = toc(tw);

            % send portions of g2 to the prior/primal nodes
            for q = 1:Qp.Value
                labSend(g2(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :), q);
            end;  g2 = []; % mem reasons
        end
    end

    %% Relative change of objective function
    % retrieve rel_x_q, norm_x_q for the workers
    rel_x = 0;
    norm_x = 0;
    for q = 1:Q
        rel_x = rel_x + rel_x_q{q};
        norm_x = norm_x + norm_x_q{q};
    end
    rel_val(t) = sqrt(rel_x / norm_x);
    end_iter(t) = toc(start_iter);

    % compute average update time (data and facet processes)
    t_facet(t) = 0; % just in case
    for q = 1:Q
        t_facet(t) = t_facet(t) + t_op{q};
    end
    t_facet(t) = t_facet(t) / Q;

    t_data(t) = 0; % just in case
    for i = Q + 1:Q + K
        t_data(t) = t_data(t) + t_op{i};
    end
    t_data(t) = t_data(t) / K;

    %% Retrieve value of the monitoring variables (residual norms + epsilons)
    norm_epsilon_check = 0;
    norm_residual_check = 0;
    for i = Q + 1:Q + K
        norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{i};
        norm_residual_check = norm_residual_check + norm_residual_check_i{i};
    end
    norm_epsilon_check = sqrt(norm_epsilon_check);
    norm_residual_check = sqrt(norm_residual_check);
    %      fprintf('Iter = %i, time = %e, t_facet = %e, t_data = %e, rel_var = %e,  epsilon = %e, residual = %e\n', t, end_iter(t), t_facet(t), t_data(t), rel_val(t), norm_epsilon_check, norm_residual_check);

    fprintf('Iter = %i, rel_var = %e,  epsilon = %e, residual = %e\n', t, rel_val(t), norm_epsilon_check, norm_residual_check);

    %% Display
    if t == 1
        fprintf('Iter = %i, time = %e, t_facet = %e, t_data = %e, rel_var = %e,  epsilon = %e, residual = %e\n', t, end_iter(t), t_facet(t), t_data(t), rel_val(t), norm_epsilon_check, norm_residual_check);
    elseif ~mod(t, 100)
        fprintf('Iter = %i, time = %e, t_facet = %e, t_data = %e, rel_var = %e,  epsilon = %e, residual = %e\n', t, end_iter(t), t_facet(t), t_data(t), rel_val(t), norm_epsilon_check, norm_residual_check);
        %% compute value of the priors in parallel
        % TODO: move this block in l.619 if computing the norsm at each iteration
        spmd
            if labindex <= Qp.Value
                % compute values for the prior terms
                % x_facet = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
                x_facet(overlap(1) + 1:end, overlap(2) + 1:end, :) = xsol_q;
                x_facet = comm2d_update_borders(x_facet, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
                [l21_norm, nuclear_norm] = fhs_compute_facet_prior(x_facet, Iq, ...
                    offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
                    offsetLq, offsetRq, crop_sparsity, crop_low_rank, apodization_window, size(v1_));
            end
        end

        % retrieve value of the priors
        l21 = 0;
        nuclear = 0;
        for q = 1:Q
            l21 = l21 + l21_norm{q};
            nuclear = nuclear + nuclear_norm{q};
        end
        % previous_obj = obj;
        % obj = (param.gamma*l21 + param.gamma0*nuclear);
        % rel_obj = abs(previous_obj - obj)/previous_obj;

        % Log
        if param.verbose >= 1
            fprintf('Iter %i\n', t);
            fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
            fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);

            if flag_synth_data
                % get xsol back from the workers
                for q = 1:Q
                    xsol(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :) = xsol_q{q};
                end
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

    %% check conditions
    % convergence pdfb (inner solver): condition
    pdfb_converged = (t - redefine_min_task_last_step_iter >= param.pdfb_min_iter) && ... % minimum number of pdfb iterations
        (t - redefine_min_task_last_step_iter >= param.pdfb_max_iter || ... % maximum number of pdfb iterations reached
        (rel_val(t) <= param.pdfb_rel_var && norm_residual_check <= param.pdfb_fidelity_tolerance * norm_epsilon_check) || ... % relative variation solution, objective and data fidelity within tolerance
        rel_val(t) <= param.pdfb_rel_var_low ... % relative variation really small and data fidelity criterion not satisfied yet
        );

    %  Update epsilons and regularisation param: condition
    pdfb_adjust_noise = (adjust_flag_noise  && ...% flag to activate adjustment of the noise level
        (norm_residual_check > param.pdfb_fidelity_tolerance * norm_epsilon_check)  && ...
        (t - redefine_min_task_last_step_iter >= param.pdfb_max_iter  || ... % trigger adjustement if  max num of pdfb iterations reached
        (redefine_min_task_last_step_iter < 1) || ... % trigger adjustement based on the itr num
        (rel_val(t) <= adjust_noise_rel_var && t - redefine_min_task_last_step_iter >= adjust_noise_min_iter)) ... % trigger adjustement if relative var. reached
        );

    if pdfb_adjust_noise || pdfb_converged
        %% update xsol
        for q = 1:Q
            xsol(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :) = xsol_q{q};
        end
        xsol_node = Composite();

        for iLab = Q + 1:numel(xsol_q)
            xsol_node{iLab} = xsol(:, :, spectral_chunk{iLab - Q});
        end
    end

    if pdfb_adjust_noise
        pdfb_converged = false;
        %% adjust noise level and resulting params.
        spmd
            if labindex > Qp.Value
                for i = 1:numel(epsilon)
                    noise_full_vect = [];
                    for j = 1:numel(epsilon{i})
                        if epsilon{i}{j} < ((2 - param.pdfb_fidelity_tolerance) * norm_res{i}{j})
                            epsilon{i}{j} = (epsilon{i}{j} + norm_res{i}{j}) * adjust_noise_change_percentage;
                        end
                        nmeas_blk = numel(y{i}{j});
                        sigma_noise_blk = sqrt(epsilon{i}{j}.^2 / nmeas_blk);
                        noise_full_vect = [noise_full_vect; ...
                            sigma_noise_blk * (randn(nmeas_blk, 1) + 1i * randn(nmeas_blk, 1)) / sqrt(2)];
                    end
                    global_sigma_noise_cmpst(i, 1) = std(noise_full_vect);
                    noise_full_vect = [];
                end
            end
        end

        % Operator's norm
        squared_operator_norm = param.squared_operator_norm(:);
        global_sigma_noise  = [];
        for i = 1:K
            global_sigma_noise =  [global_sigma_noise; global_sigma_noise_cmpst{Q + i}];
        end
        % noise level / regularization parameter
        [sig, sig_bar, mu_chi, sig_chi, sig_sara] = ...
            compute_noise_level(M, N, K, global_sigma_noise(:), ...
            'fhs', Qx, Qy, overlap_size, squared_operator_norm);
        % apply multiplicative factor for the regularization parameters (if needed)
        gamma0 = alph_bar * sig_bar;
        gamma = alph * sig;
        beta0 = Composite();
        for q = 1:Q
            beta0{q} = gamma0(q) / sigma0;
        end
        beta1 = gamma / sigma1; %         beta1 = parallel.pool.Constant(gamma / sigma1);

        redefine_min_task_last_step_iter = t;

        fprintf('mu_chi = %.4e, sig_chi = %.4e, sig_sara = %.4e\n', mu_chi, sig_chi, sig_sara);
        fprintf('Noise levels: sig = %.4e, sig_bar = [%.4e, %.4e]\n', sig, min(sig_bar), max(sig_bar));
        fprintf('Additional multiplicative actors gam = %.4e, gam_bar = %.4e\n', alph, alph_bar);
        fprintf('Regularization parameters: mu = %.4e, mu_bar = %.4e\n', gamma, mean(gamma0));
        try  fitswrite(xsol, [checkpoint_name '_tmp_wb_model_noise_adjust.fits']);
        end
    elseif pdfb_converged
        adjust_flag_noise = false;
        %% Reweighting (in parallel)
        % Evaluate relative variation for the reweighting scheme
        spmd
            if labindex <= Qp.Value
                rel_x_reweighting_q = norm(xlast_reweight_q(:) - xsol_q(:))^2;
                norm_x_reweighting_q = norm(xlast_reweight_q(:))^2;
                xlast_reweight_q = []; % AD: for mem reasons, will be updated later.
                % xlast_reweight_q = xsol_q;
            end
        end
        rel_x_reweighting = 0;
        norm_x_reweighting = 0;
        for q = 1:Q
            rel_x_reweighting = rel_x_reweighting + rel_x_reweighting_q{q};
            norm_x_reweighting = norm_x_reweighting + norm_x_reweighting_q{q};
        end
        rel_x_reweighting = sqrt(rel_x_reweighting / norm_x_reweighting);

        reweighting_converged =  ~reweighting_flag || ...
            (pdfb_converged && ...                  % do not exit solver before the current pdfb algorithm converged
            reweight_step_count >= param.reweighting_min_iter && ...   % minimum number of reweighting iterations
            (reweight_step_count >= param.reweighting_max_iter || ... % maximum number of reweighting iterations reached
            rel_x_reweighting <= param.reweighting_rel_var) ...         % relative variation
            );

        if reweighting_converged
            flag_convergence = 1;
            break
        end

        fprintf('Reweighting: %i, relative variation: %e, reweighting parameter: %e \n\n', reweight_step_count + 1, rel_x_reweighting, reweighting_alpha);

        spmd
            if labindex <= Qp.Value
                % update weights
                fprintf('\ngetting facets');
                x_facet(overlap(1) + 1:end, overlap(2) + 1:end, :) = xsol_q;
                x_facet = comm2d_update_borders(x_facet, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);

                fprintf('\nUpdating weights,  xfacet size  %d x  %d', size(x_facet, 1), size(x_facet, 2));
                % ! -- New reweighting with proper floor level
                [weights1_, weights0_] = fhs_update_weights(x_facet, size(v1_), ...
                    Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, ...
                    Ncoefs_q, dims_overlap_ref_q, offsetLq, offsetRq, ...
                    reweighting_alphap, crop_sparsity, crop_low_rank, apodization_window, sig_, sig_bar_);
                % ! HOMOTOPY (deactivated)
                % if flag_homotopy
                %     reweighting_alphap = max(reweighting_alpha_ffp.Value * reweighting_alphap, 1);
                % end
                % ! --
            else
                % compute residual image on the data nodes
                % res_ = compute_residual_images(xsol(:, :, spectral_chunk{labindex - Qp.Value}), y, A, At, G, W, flag_dimensionality_reduction, Lambda);
                res_ = compute_residual_images(xsol_node, y, A, At, G, W, flag_dimensionality_reduction, Lambda);
                xsol_node = []; % mem reasons
                fprintf('\nResidual computed');
            end
        end; clear xsol_node;
        % ! -- HOMOTOPY (deactivated)
        % if flag_homotopy
        %     reweighting_alpha = max(param.reweighting_alpha_ff * reweighting_alpha, 1);
        % end
        % ! --
        param.reweighting_alpha = reweighting_alpha;
        param.init_reweight_step_count = reweight_step_count + 1;
        param.init_reweight_last_iter_step = t;
        param.init_t_start = t + 1;

        %% compute value of the priors in parallel
        spmd
            if labindex <= Qp.Value
                % compute values for the prior terms
                % x_facet = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
                x_facet(overlap(1) + 1:end, overlap(2) + 1:end, :) = xsol_q;
                x_facet = comm2d_update_borders(x_facet, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
                [l21_norm, nuclear_norm] = fhs_compute_facet_prior(x_facet, Iq, ...
                    offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
                    offsetLq, offsetRq, crop_sparsity, crop_low_rank, apodization_window, size(v1_));
            end
        end

        % keep record of last sol
        spmd
            if labindex <= Qp.Value
                xlast_reweight_q = xsol_q;
            end
        end

        % retrieve value of the priors
        l21 = 0;
        nuclear = 0;
        for q = 1:Q
            l21 = l21 + l21_norm{q};
            nuclear = nuclear + nuclear_norm{q};
        end
        % obj = param.gamma0*nuclear + param.gamma*l21;

        if flag_synth_data
            % get xsol back from the workers

            sol = reshape(xsol(:), numel(xsol(:)) / n_channels, n_channels);
            SNR = 20 * log10(norm(X0(:)) / norm(X0(:) - sol(:)));
            psnrh = zeros(n_channels, 1);
            for i = 1:n_channels
                psnrh(i) = 20 * log10(norm(X0(:, i)) / norm(X0(:, i) - sol(:, i)));
            end
            SNR_average = mean(psnrh);
            fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
        end

        if (reweight_step_count == 0) || (reweight_step_count == 1) || ...
                (~mod(reweight_step_count, param.backup_frequency))
            % save model estimate
            fitswrite(xsol, [checkpoint_name '_CURRENT_WB_MODEL.fits']);
            for k = 1:K;  res(:, :, spectral_chunk{k}) = res_{Q + k};
            end
            fitswrite(res, [checkpoint_name '_CURRENT_WB_RESIDUAL.fits']);
            clear res;
            if param.save_intermediate_results_mat
                m = matfile(strcat(checkpoint_name, '_rw', num2str(reweighting_flag), '.mat'), ...
                    'Writable', true);
                m.param = param;
                m.res = zeros(size(xsol));
                m.g = zeros(size(xsol));
                m.xsol = zeros(size(xsol));
                m.epsilon = cell(K, 1);
                m.v2 = cell(K, 1);
                m.proj = cell(K, 1);
                m.t_block = cell(K, 1);
                m.norm_res = cell(K, 1);
                m.v0 = cell(Q, 1);
                m.v1 = cell(Q, 1);
                m.weights0 = cell(Q, 1);
                m.weights1 = cell(Q, 1);
                % Retrieve variables from workers
                % facet nodes
                for q = 1:Q
                    m.v0(q, 1) = v0_(q);
                    m.v1(q, 1) = v1_(q);
                    m.weights0(q, 1) = weights0_(q);
                    m.weights1(q, 1) = weights1_(q);
                    %                 m.xsol(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :) = xsol_q{q};
                    m.g(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :) = g_q{q};
                end
                m.xsol = xsol; % xsol_q{q};

                % data nodes
                for k = 1:K
                    m.res(:, :, spectral_chunk{k}) = res_{Q + k};
                    res_{Q + k} = [];
                    m.v2(k, 1) = v2_(Q + k);
                    m.proj(k, 1) = proj_(Q + k);
                    m.t_block(k, 1) = t_block(Q + k);
                    m.epsilon(k, 1) = epsilon(Q + k);
                    m.norm_res(k, 1) = norm_res(Q + k);
                end
                m.end_iter = end_iter;
                m.t_facet = t_facet;
                m.t_data = t_data;
                m.rel_val = rel_val;

                if flag_synth_data
                    m.SNR = SNR;
                    m.SNR_average = SNR_average;
                end
                clear m;
            end
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
        redefine_min_task_last_step_iter = t;
        if reweight_step_count >= param.reweighting_max_iter
            fprintf('\n\n No more reweights \n\n');
        end
    end
end
toc(start_loop);

% Collect image facets back to the master
for q = 1:Q
    xsol(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :) = xsol_q{q};
end

% Calculate residual images
xsol_node = Composite();
for iLab = Q + 1:numel(xsol_q)
    xsol_node{iLab} = xsol(:, :, spectral_chunk{iLab - Q});
end
spmd
    if labindex <= Qp.Value
        x_facet(overlap(1) + 1:end, overlap(2) + 1:end, :) = xsol_q;
        x_facet = comm2d_update_borders(x_facet, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);

        [l21_norm, nuclear_norm] = fhs_compute_facet_prior(x_facet, Iq, ...
            offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
            offsetLq, offsetRq, crop_sparsity, crop_low_rank, apodization_window, size(v1_));
    else
        % res_ = compute_residual_images(xsol(:, :, spectral_chunk{labindex - Qp.Value}), y, A, At, G, W, flag_dimensionality_reduction, Lambda);
        res_ = compute_residual_images(xsol_node, y, A, At, G, W, flag_dimensionality_reduction, Lambda); xsol_node = []; % mem reasons
    end
end; clear xsol_node;
%  fits
m = matfile([checkpoint_name, '_rw' num2str(reweighting_flag) '.mat'], ...
    'Writable', true);
m.param = param;
m.res = zeros(size(xsol));
m.g = zeros(size(xsol));
m.xsol = zeros(size(xsol));
m.epsilon = cell(K, 1);
m.v2 = cell(K, 1);
m.proj = cell(K, 1);
m.t_block = cell(K, 1);
m.norm_res = cell(K, 1);
m.v0 = cell(Q, 1);
m.v1 = cell(Q, 1);
m.weights0 = cell(Q, 1);
m.weights1 = cell(Q, 1);

% Retrieve variables from workers
% facet nodes
for q = 1:Q
    m.v0(q, 1) = v0_(q);
    v0_{q} = [];
    m.v1(q, 1) = v1_(q);
    v1_{q} = [];
    m.weights0(q, 1) = weights0_(q);
    weights0_{q} = [];
    m.weights1(q, 1) = weights1_(q);
    weights1_{q} = [];
    m.g(I(q, 1) + 1:I(q, 1) + dims(q, 1), I(q, 2) + 1:I(q, 2) + dims(q, 2), :) = g_q{q};
    g_q{q} = [];
end

% data nodes
for k = 1:K
    m.res(:, :, spectral_chunk{k}) = res_{Q + k};
    res_{Q + k} = [];
    m.v2(k, 1) = v2_(Q + k);
    v2_{Q + k} = [];
    m.proj(k, 1) = proj_(Q + k);
    proj_{Q + k} = [];
    m.t_block(k, 1) = t_block(Q + k);
    t_block{Q + k} = [];
    m.epsilon(k, 1) = epsilon(Q + k);
    epsilon{Q + k} = [];
    m.norm_res(k, 1) = norm_res(Q + k);
end
m.xsol = xsol;
% norm_res_out = sqrt(sum(sum(sum((m.res).^2))));

% Update param structure and save
param.reweighting_alpha = reweighting_alpha;
param.init_reweight_step_count = reweight_step_count;
param.init_reweight_last_iter_step = t;
param.init_t_start = t + 1;
m.param = param;
m.end_iter = end_iter;
m.t_facet = t_facet;
m.t_data = t_data;
m.rel_val = rel_val;
fitswrite(m.xsol, [checkpoint_name '_FINAL_WB_MODEL' '.fits']);
fitswrite(m.res, [checkpoint_name '_FINAL_WB_RESIDUAL' '.fits']);
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
l21 = 0;
nuclear = 0;
for q = 1:Q
    l21 = l21 + l21_norm{q};
    nuclear = nuclear + nuclear_norm{q};
end

norm_epsilon_check = 0;
norm_residual_check = 0;
for k = Q + 1:Q + K
    norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{k};
    norm_residual_check = norm_residual_check + norm_residual_check_i{k};
end
norm_epsilon_check = sqrt(norm_epsilon_check);
norm_residual_check = sqrt(norm_residual_check);

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
