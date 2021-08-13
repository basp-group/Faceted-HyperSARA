function xsol = facetHyperSARA(y, epsilon, A, At, pU, G, W, param, ...
    Qx, Qy, K, wavelet, filter_length, nlevel, window_type, ...
    spectral_chunk, nChannels, overlap_size, alph, alph_bar, ...
    M, N, oy, ox, name_warmstart, name_checkpoint, flagDR, Sigma, ...
    varargin)

% TODO: try to replace A, At, pU, G, W by a functor (if possible)

%facetHyperSARA_cw: faceted HyperSARA
%
% version with a fixed overlap for the faceted nuclear norm, larger or 
% smaller than the extension needed for the 2D segmented discrete wavelet 
% transforms (sdwt2). Includes spatial weihting correction for the faceted
% nuclear norm (triangular, hamming, piecewise_constant, no weights by
% default).
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
% > Qx          number of facets along dimension x [1]
% > Qy          number of facets along dimension y [1]
% > K           number of Matlab data fidelity processes [1]
% > wavelet     wavelet doctionaries considered (should contain 'self' by
%               default in last position)
% > filter_length  size of the wavelet filters considered (by cinvention, 0 for the Dirac basis)
% > nlevel      decomposition depth [1]
% > spectral_chunk    indices of the bands handled by each data node {K, 1}
% > nChannels           total number of spectral channels [1]
% > d           size of the fixed overlap for the faceted nuclear norm 
%               [1, 2]
% > window_type type of apodization window affecting the faceted nuclear
%               norm prior [string]
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
% term in a single place. Constant overlap for the nuclear norm assuming d
% is smaller than the smallest overlap for the sdwt2 (the other option 
% would also change the communication process (borders and reduction 
% operation)). d <= (power(2, nlevel)-1)*(max(filter_length(:)-1))

%% NOTE:
% this version relies on a specialised version of sdwt2, slightly less
% general but faster (based on Arwa's work).

% This function solves:
%
% min || X ||_* + lambda * ||Psit(X)||_2,1   s.t.  || Y - A(X) ||_2 <= epsilon and x>=0
%
%%

% size of the oversampled Fourier space (vectorized)
No = M*oy*N*ox;

% -- instantiate auxiliary variables for sdwt2
% define reference 2D facets (no overlap)
Q = Qx*Qy;
rg_y = split_range(Qy, M);
rg_x = split_range(Qx, N);
I = zeros(Q, 2);
dims = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        I(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1];
        dims(q, :) = [rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
    end
end
clear rg_y rg_x;

rg_yo = split_range(Qy, M, overlap_size(1));
rg_xo = split_range(Qx, N, overlap_size(2));
Io = zeros(Q, 2);
dims_o = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        Io(q, :) = [rg_yo(qy, 1)-1, rg_xo(qx, 1)-1];
        dims_o(q, :) = [rg_yo(qy,2)-rg_yo(qy,1)+1, rg_xo(qx,2)-rg_xo(qx,1)+1];
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
c_chunksp = parallel.pool.Constant(spectral_chunk);
waveletp = parallel.pool.Constant(wavelet);
nlevelp = parallel.pool.Constant(nlevel);
offsetp = parallel.pool.Constant(offset);

% define auxiliary composite variables (local to a given worker)
[Iq, dims_q, dims_oq, dims_overlap_ref_q, I_overlap_q, ...
    dims_overlap_q, status_q, offsetLq, offsetRq, Ncoefs_q, temLIdxs_q, ...
    temRIdxs_q, overlap_g_south, overlap_g_east, overlap_g_south_east, overlap, ...
    w, crop_nuclear, crop_l21] = setup_priors(Qx, Qy, I, dims, dims_o, ...
    dims_overlap_ref, I_overlap, dims_overlap, status, offsetL, offsetR, ...
    Ncoefs, temLIdxs, temRIdxs, window_type, overlap_size);

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
%! --

g_q = Composite();
xsol_q = Composite();
if init_flag
   size(xsol)
       	for q = 1:Q
        xsol_q{q} = xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :);
        g_q{q} = init_m.g(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :);
    end
    fprintf('g uploaded \n\n')
else
    for q = 1:Q
        xsol_q{q} = xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :);
        g_q{q} = zeros([dims(q, :), nChannels]);
    end
    fprintf('g initialized \n\n')
end

%! -- TO BE CHECKED
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
reweighting_alpha_ffp = parallel.pool.Constant(param.reweighting_alpha_ff);

if isfield(param,'init_reweight_step_count')
    reweight_step_count = param.init_reweight_step_count;
    fprintf('reweight_step_count uploaded\n\n')
else
    param.init_reweight_step_count = 0;
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
% ! assumes backup file exactly saved at the time a reweighting step occured
% ! (initialized to xsol_q)
v0_ = Composite();
weights0_ = Composite();
v1_ = Composite();
weights1_ = Composite();
xlast_reweight_q = Composite(); 
if init_flag
    for q = 1:Q
        v0_(q) = init_m.v0(q,1);
        v1_(q) = init_m.v1(q,1);
        weights0_(q) = init_m.weights0(q,1);
        weights1_(q) = init_m.weights1(q,1);
    end
    spmd
        if labindex <= Qp.Value
            max_dims = max(dims_overlap_ref_q, dims_oq);
            xlast_reweight_q = xsol_q;
        end
    end
    fprintf('v0, v1, weigths0, weights1 uploaded \n\n')
else
    spmd
        if labindex <= Qp.Value
            max_dims = max(dims_overlap_ref_q, dims_oq);
            xlast_reweight_q = xsol_q;
            %!-- TO BE CHECKED           
            x_overlap = zeros([max_dims, size(xsol_q, 3)]);
            x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
            x_overlap = comm2d_update_borders(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);

            % weights initialized from initial primal variable, dual variables to 0
            [v0_, v1_, weights0_, weights1_] = initialize_dual_and_weights(x_overlap, ...
                Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, max_dims-crop_nuclear, nChannels, dims_overlap_ref_q, ...
                offsetLq, offsetRq, reweighting_alphap, crop_l21, crop_nuclear, w, sig_, sig_bar_);

            % weights and dual variables initialized from initial primal variable
            % [v0, v1, weights0, weights1] = initialize_dual_and_weights2(x_overlap, ...
            %     Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
            %     offsetLq, offsetRq, reweighting_alphap, crop_l21, crop_nuclear, w);
            %!--
        end
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
adapt_eps_change_percentage = parallel.pool.Constant(param.adapt_eps_change_percentage);

if init_flag
    norm_res = Composite();
    v2_ = Composite();
    t_block = Composite();
    proj_ = Composite();
    for k = 1:K
        norm_res(Q+k) = init_m.norm_res(k,1);
        v2_(Q+k) = init_m.v2(k,1);
        proj_(Q+k) = init_m.proj(k,1);
        t_block(Q+k) = init_m.t_block(k,1);
    end
    fprintf('v2, proj, t_block, norm_res uploaded \n\n')
else
    % ! assumes primal variable initialized to 0
    spmd
        if labindex > Qp.Value
            [v2_, norm_res, t_block, proj_] = initilize_data_worker(y);
        end
    end
    fprintf('v2, proj, t_block, norm_res initialized \n\n')
end

% initialize xi and Fxi_old
spmd
    if labindex <= Qp.Value
        % send xhat_q (communication towards the data nodes)
        for i = 1:K
            labSend(xsol_q(:,:,c_chunksp.Value{i}), Qp.Value+i);
        end
    else
        xi = zeros(M, N, numel(y));
        for q = 1:Qp.Value
            xi(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2), :) = ...
                labReceive(q);
        end

        Fxi_old = zeros(No, numel(y));
        for l = 1:numel(y)
            Fxi_old(:, l) = A(xi(:,:,l));
        end
    end
end

%Step sizes for the dual variables
sigma0 = 1.0/param.nu0;
sigma1 = 1.0/param.nu1;
sigma2 = 1.0./param.nu2;

% Step size primal
tau = 0.99/3;

% Update constant dual variables
sigma00 = parallel.pool.Constant(tau*sigma0);
sigma11 = parallel.pool.Constant(tau*sigma1);
% sigma22 = parallel.pool.Constant(tau*sigma2);
sigma22 = Composite();
for k = 1:K
    sigma22{Q+k} = tau*sigma2(spectral_chunk{k});
end
% beta0 = parallel.pool.Constant(param.gamma0/sigma0);
beta0 = Composite();
for q = 1:Q
    beta0{q} = param.gamma0(q)/sigma0;
end
beta1 = parallel.pool.Constant(param.gamma/sigma1);
param.alph = alph;
param.alph_bar = alph_bar;

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

max_iter = (param.reweighting_max_iter + 1)*param.pdfb_max_iter;
if init_flag
    rel_val = init_m.rel_val;
    end_iter = init_m.end_iter;
    t_facet = init_m.t_facet;
    t_data = init_m.t_data;
    fprintf('rel_val, end_iter, t_facet and t_data uploaded \n\n')
else
    rel_val = zeros(max_iter, 1);
    end_iter = zeros(max_iter, 1);
    t_facet = zeros(max_iter, 1);
    t_data = zeros(max_iter, 1);
    fprintf('rel_val, end_iter, t_facet and t_data initialized \n\n')
end

%! check warm-start worked as expected
if init_flag
    spmd
        if labindex > Qp.Value
            [norm_residual_check_i, norm_epsilon_check_i] = sanity_check(epsilon, norm_res);
        end
    end
    norm_epsilon_check = 0;
    norm_residual_check = 0;
    for i = Q+1:Q+K
        norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{i};
        norm_residual_check = norm_residual_check + norm_residual_check_i{i};
    end
    norm_epsilon_check = sqrt(norm_epsilon_check);
    norm_residual_check = sqrt(norm_residual_check);

    % compute value of the priors in parallel
    spmd
        if labindex <= Qp.Value
            % compute values for the prior terms
            x_overlap = zeros([max_dims, size(xsol_q, 3)]);
            x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
            x_overlap = comm2d_update_borders(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
            [l21_norm, nuclear_norm] = compute_facet_prior_overlap(x_overlap, Iq, ...
                offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
                offsetLq, offsetRq, crop_l21, crop_nuclear, w, size(v1_));
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
    if (param.verbose >= 1)
        fprintf('Iter %i\n',max(t_start-1, 1));
        fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(max(t_start-1, 1)));
        fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);

        if flag_synth_data
            for q = 1:Q
                xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
            end
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

start_loop = tic;

fprintf('START THE LOOP MNRAS ver \n\n')

for t = t_start : max_iter
    
    start_iter = tic;
    
    spmd
        if labindex <= Qp.Value
            % primal/prior nodes (1:Q)
            
            % update primal variable
            tw = tic;
            [xsol_q, xhat_q, rel_x_q, norm_x_q] = update_primal(xsol_q, g_q);
            t_op = toc(tw);
            
            % send xhat_q (communication towards the data nodes)
            for k = 1:K
                labSend(xsol_q(:,:,c_chunksp.Value{k}), Qp.Value+k);
            end
            
            % update borders (-> versions of xhat with overlap)
            tw = tic;
            x_overlap = zeros([max_dims, size(xsol_q, 3)]);
            x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xhat_q;
            x_overlap = comm2d_update_borders(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
            
            % update dual variables (nuclear, l21)
            [v0_, g0] = update_dual_nuclear(v0_, x_overlap(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :), w, weights0_, beta0);
            [v1_, g1] = update_dual_l21(v1_, x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :), weights1_, beta1.Value, Iq, ...
                dims_q, I_overlap_q, dims_overlap_q, offsetp.Value, status_q, ...
                nlevelp.Value, waveletp.Value, Ncoefs_q, temLIdxs_q, temRIdxs_q, offsetLq, offsetRq, dims_overlap_ref_q);
            g = zeros(size(x_overlap));
            g(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :) = sigma00.Value*g0;
            g(crop_l21(1)+1:end, crop_l21(2)+1:end, :) = g(crop_l21(1)+1:end, crop_l21(2)+1:end, :) + sigma11.Value*g1;          
            g = comm2d_reduce(g, overlap, Qyp.Value, Qxp.Value);
            t_op = t_op + toc(tw);
            
            % compute g_ for the final update term
            g_q = g(overlap(1)+1:end, overlap(2)+1:end, :);
            
            % retrieve portions of g2 from the data nodes
            for k = 1:Kp.Value
                g_q(:,:,c_chunksp.Value{k}) = g_q(:,:,c_chunksp.Value{k}) + labReceive(Qp.Value+k);
            end
        else
            % data nodes (Q+1:Q+K) (no parallelisation over data blocks, just frequency)
            % retrieve xhat_i from the prior/primal nodes
            for q = 1:Qp.Value
                xi(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2), :) = ...
                    labReceive(q);
            end
            tw = tic;
            [v2_, g2, Fxi_old, proj_, norm_res, norm_residual_check_i, norm_epsilon_check_i] = update_dual_fidelity(v2_, y, xi, ...
                Fxi_old, proj_, A, At, G, W, pU, epsilon, ...
                elipse_proj_max_iter.Value, elipse_proj_min_iter.Value, ...
                elipse_proj_eps.Value, sigma22, flagDR, Sigma);
            t_op = toc(tw);

            % send portions of g2 to the prior/primal nodes
            for q = 1:Qp.Value
                labSend(g2(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2), :), q);
            end
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
    rel_val(t) = sqrt(rel_x/norm_x);
    end_iter(t) = toc(start_iter);

    % compute average update time (data and facet processes)
    t_facet(t) = 0; % just in case
    for q = 1:Q
        t_facet(t) = t_facet(t) + t_op{q};
    end
    t_facet(t) = t_facet(t)/Q;

    t_data(t) = 0; % just in case
    for i = Q+1:Q+K
        t_data(t) = t_data(t) + t_op{i};
    end
    t_data(t) = t_data(t)/K;
    
    %% Retrieve value of the monitoring variables (residual norms + epsilons)
    norm_epsilon_check = 0;
    norm_residual_check = 0;
    for i = Q+1:Q+K
        norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{i};
        norm_residual_check = norm_residual_check + norm_residual_check_i{i};
    end
    norm_epsilon_check = sqrt(norm_epsilon_check);
    norm_residual_check = sqrt(norm_residual_check);

    fprintf('Iter = %i, time = %e, t_facet = %e, t_data = %e, rel_var = %e,  epsilon = %e, residual = %e\n',t,end_iter(t),t_facet(t),t_data(t),rel_val(t),norm_epsilon_check,norm_residual_check);
    
    %% Display
    if ~mod(t,100)
        
        %% compute value of the priors in parallel
        % TODO: move this block in l.619 if computing the norsm at each iteration
        spmd
            if labindex <= Qp.Value
                % compute values for the prior terms
                %x_overlap = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
                x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
                x_overlap = comm2d_update_borders(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
                [l21_norm, nuclear_norm] = compute_facet_prior_overlap(x_overlap, Iq, ...
                    offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
                    offsetLq, offsetRq, crop_l21, crop_nuclear, w, size(v1_));
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
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
            fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);

            if flag_synth_data
                % get xsol back from the workers
                for q = 1:Q
                    xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
                end
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
    % && rel_obj <= param.pdfb_rel_obj

    %% Update epsilons (in parallel)
    flag_epsilonUpdate = param.use_adapt_eps && ...  % activate espilon update 
    (t > param.adapt_eps_start) && ...               % update allowed after a minimum of iterations in the 1st reweighting
    (rel_val(t) < param.adapt_eps_rel_var);          % relative variation between 2 consecutive pdfb iterations
    
    if flag_epsilonUpdate
        spmd
            if labindex > Qp.Value
                [epsilon, t_block] = update_epsilon(epsilon, t, t_block, norm_res, ...
                    adapt_eps_tol_in.Value, adapt_eps_tol_out.Value, adapt_eps_steps.Value, ...
                    adapt_eps_change_percentage.Value);
            end
        end
    end
    %! --
    
    %% Reweighting (in parallel)
    if pdfb_converged
        % Evaluate relative variation for the reweighting scheme
        spmd 
            if labindex <= Qp.Value
                rel_x_reweighting_q = norm(xlast_reweight_q(:) - xsol_q(:))^2;
                norm_x_reweighting_q = norm(xlast_reweight_q(:))^2;
                xlast_reweight_q = xsol_q;
            end
        end
        rel_x_reweighting = 0;
        norm_x_reweighting = 0;    
        for q = 1:Q
            rel_x_reweighting = rel_x_reweighting + rel_x_reweighting_q{q};
            norm_x_reweighting = norm_x_reweighting + norm_x_reweighting_q{q};
        end
        rel_x_reweighting = sqrt(rel_x_reweighting/norm_x_reweighting);

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
        
        spmd
            if labindex <= Qp.Value
                % update weights
                x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
                x_overlap = comm2d_update_borders(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);

                % [weights1_, weights0_] = update_weights_overlap(x_overlap, size(v1_), ...
                %     Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, ...
                %     Ncoefs_q, dims_overlap_ref_q, offsetLq, offsetRq, ...
                %     reweighting_alphap, crop_l21, crop_nuclear, w);

                %! -- TO BE CHECKED (using new reweighting with proper floor level)
                [weights1_, weights0_] = update_weights_overlap2(x_overlap, size(v1_), ...
                Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, ...
                Ncoefs_q, dims_overlap_ref_q, offsetLq, offsetRq, ...
                reweighting_alphap, crop_l21, crop_nuclear, w, sig_, sig_bar_);
                if flag_homotopy
                    reweighting_alphap = max(reweighting_alpha_ffp.Value*reweighting_alphap, 1);
                end
                %! --
            else
                % compute residual image on the data nodes
                res_ = compute_residual_images(xsol(:,:,spectral_chunk{labindex-Qp.Value}), y, A, At, G, W, flagDR, Sigma);
            end
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

        %% compute value of the priors in parallel
        spmd
            if labindex <= Qp.Value
            % compute values for the prior terms
            %x_overlap = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
            x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
            x_overlap = comm2d_update_borders(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
            [l21_norm, nuclear_norm] = compute_facet_prior_overlap(x_overlap, Iq, ...
                offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
                offsetLq, offsetRq, crop_l21, crop_nuclear, w, size(v1_));
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
            for q = 1:Q
                xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
            end
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
                m.v0(q,1) = v0_(q);
                m.v1(q,1) = v1_(q);
                m.weights0(q,1) = weights0_(q);
                m.weights1(q,1) = weights1_(q);
                m.xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
                m.g(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = g_q{q};
            end
            % data nodes
            for k = 1:K
                m.res(:,:,spectral_chunk{k}) = res_{Q+k};
                res_{Q+k} = [];
                m.v2(k,1) = v2_(Q+k);
                m.proj(k,1) = proj_(Q+k);
                m.t_block(k,1) = t_block(Q+k);
                m.epsilon(k,1) = epsilon(Q+k);
                m.norm_res(k,1) = norm_res(Q+k);
            end
            m.end_iter = end_iter;
            m.t_facet = t_facet;
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
                if flag_synth_data
                    fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
                end
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

% Collect image facets back to the master
for q = 1:Q
    xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
end

% Calculate residual images
spmd
    if labindex <= Qp.Value
        x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
        x_overlap = comm2d_update_borders(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
        
        [l21_norm, nuclear_norm] = compute_facet_prior_overlap(x_overlap, Iq, ...
            offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
            offsetLq, offsetRq, crop_l21, crop_nuclear, w, size(v1_));
    else
        res_ = compute_residual_images(xsol(:,:,spectral_chunk{labindex-Qp.Value}), y, A, At, G, W, flagDR, Sigma);
    end
end

m = matfile([name_checkpoint, '_rw=' num2str(reweight_step_count) '.mat'], ...
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
    m.v0(q,1) = v0_(q);
    v0_{q} = [];
    m.v1(q,1) = v1_(q);
    v1_{q} = [];
    m.weights0(q,1) = weights0_(q);
    weights0_{q} = [];
    m.weights1(q,1) = weights1_(q);
    weights1_{q} = [];
    m.g(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = g_q{q};
    g_q{q} = [];
end

% data nodes
for k = 1:K
    m.res(:,:,spectral_chunk{k}) = res_{Q+k};
    res_{Q+k} = [];
    m.v2(k,1) = v2_(Q+k);
    v2_{Q+k} = [];
    m.proj(k,1) = proj_(Q+k);
    proj_{Q+k} = [];
    m.t_block(k,1) = t_block(Q+k);
    t_block{Q+k} = [];
    m.epsilon(k,1) = epsilon(Q+k);
    epsilon{Q+k} = [];
    m.norm_res(k,1) = norm_res(Q+k);
end
m.xsol = xsol;
% norm_res_out = sqrt(sum(sum(sum((m.res).^2))));

% Update param structure and save
param.reweighting_alpha = reweighting_alpha;
param.init_reweight_step_count = reweight_step_count;
param.init_reweight_last_iter_step = t;
param.init_t_start = t+1;
m.param = param;
m.end_iter = end_iter;
m.t_facet = t_facet;
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
l21 = 0;
nuclear = 0;
for q = 1:Q
    l21 = l21 + l21_norm{q};
    nuclear = nuclear + nuclear_norm{q};
end

norm_epsilon_check = 0;
norm_residual_check = 0;
for k = Q+1:Q+K
    norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{k};
    norm_residual_check = norm_residual_check + norm_residual_check_i{k};
end
norm_epsilon_check = sqrt(norm_epsilon_check);
norm_residual_check = sqrt(norm_residual_check);

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
