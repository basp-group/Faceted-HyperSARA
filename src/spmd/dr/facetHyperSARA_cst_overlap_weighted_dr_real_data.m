function [xsol,param,epsilon,t,rel_val,nuclear,l21,norm_res_out,end_iter] = ...
    facetHyperSARA_cst_overlap_weighted_dr_real_data(y, imsize, ...
    epsilon, A, At, H, pU, T, W, param, Qx, Qy, K, wavelet, L, nlevel, ...
    c_chunks, c, d, window_type, init_file_name, name)
%facetHyperSARA_cst_overlap_weighted_dr_real_data: faceted HyperSARA
%
% version with a fixed overlap for the faceted nuclear norm, larger or 
% smaller than the extension needed for the 2D segmented discrete wavelet 
% transforms (sdwt2). Includes spatial weihting correction for the faceted
% nuclear norm (triangular, hamming, piecewise_constant, no weights by
% default). Leverages dimensionality redution (DR). Version for real data.
%
%-------------------------------------------------------------------------%
%%
% Input: 
%
% > y           blocks of visivilities {L}{nblocks_l}
% > imsize      size of the wideband image [1, 2]
% > epsilon     l2-ball norms {L}{nblocks_l}
% > A           measurement operator
% > At          adjoint of the measurement operator
% > H           holographic matrices G'*G {L}{nblocks_l}
% > pU          preconditioning matrices {L}{nblocks_l}
% > T           pseudo singular values from the reduction operator 
%               {L}{nblocks_l}
% > W           masks for selection of the blocks of visibilities
% > param       algorithm parameters (struct)
%
%   general
%   > .verbose           print log or not
%   > .rel_var   (1e-5)  stopping criterion
%   > .max_iter (10000)  max number of iterations
%
%   convergence
%   > .nu0 = 1
%   > .nu1      upper bound on the norm of the operator Psi
%   > .nu2      norm of the faceting operator (= 1)
%   > .gamma    regularization parameter (l21 norm)
%
%   reweighting
%   > .reweight_alpha_ff   (0.9)        
%   > .total_reweights      (30)     -1 if you don't want reweighting
%   > .use_reweight_steps    (1)     reweighting by fixed steps
%   > .reweight_step_size  (300)     reweighting step size
%   > .reweight_steps = [5000: param_HSI.reweight_step_size :10000];
%   > .step_flag = 1;
%   > .use_reweight_eps   (false)    reweighting w.r.t the relative change of the solution
%   > .reweight_max_reweight_itr     param_HSI.max_iter - param_HSI.reweight_step_size;
%   > .reweight_rel_var   (1e-4)     criterion for performing reweighting
%   > .reweight_min_steps_rel_var (300) minimum number of iterations between consecutive reweights
%   > .sig                           noise level (wavelet space)
%   > .sig_bar                       noise level (singular value space)

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
%
% > Qx          number of facets along dimension x [1]
% > Qy          number of facets along dimension y [1]
% > K           number of datab computing processes [1]
% > wavelet     wavelet doctionaries considered (should contain 'self' by
%               default in last position)
% > L           size of the wavelet filters considered (by cinvention, 0 for the Dirac basis)
% > nlevel      decomposition depth [1]
% > c_chunks    indices of the bands handled by each data node {K, 1}
% > c           total number of spectral channels [1]
% > d           size of the fixed overlap for the faceted nuclear norm 
%               [1, 2]
% > window_type type of apodization window affecting the faceted nuclear
%               norm prior [string]
% > init_file_name  name of a valid .mat file for initialization (for warm-restart)
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
% < reweight_alpha  last value of the reweigthing parameter [1]
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
% Code: P.-A. Thouvenin, M. Jiang, A. Abdulaziz
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
% would also change the communication process (ghost cells and reduction 
% operation)). d <= (power(2, nlevel)-1)*(max(L(:)-1))

%% NOTE:
% this version relies on a specialised version of sdwt2, slightly less
% general but faster (based on Arwa's work).

% This function solves:
%
% min || X ||_* + lambda * ||Psit(X)||_2,1   s.t.  || Y - A(X) ||_2 <= epsilon and x>=0
%
% Author: P.-A. Thouvenin.
% Adapted from codes by: Abdullah Abdulaziz.

%% REMARKS
% ...
%%

% maxNumCompThreads(param.num_workers);

% initialize monitoring variables (display active)
norm_epsilon_check = Inf;
norm_residual_check = 0;

% size of the oversampled Fourier space (vectorized)
No = size(W{1}{1}{1}, 1);

% number of pixels (spatial dimensions)
M = imsize(1);
N = imsize(2);

% define reference spatial facets (no overlap)
Q = Qx*Qy;
rg_y = domain_decomposition(Qy, M);
rg_x = domain_decomposition(Qx, N);
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

%%- begin initialization sdwt2
% instantiate auxiliary variables for sdwt2
[~, ~, I_overlap_ref, dims_overlap_ref, I_overlap, dims_overlap, ...
    ~, ~, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices([M, N], I, dims, nlevel, wavelet, L);

rg_yo = domain_decomposition_overlap2(Qy, M, d);
rg_xo = domain_decomposition_overlap2(Qx, N, d);
Io = zeros(Q, 2);
dims_o = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        Io(q, :) = [rg_yo(qy, 1)-1, rg_xo(qx, 1)-1];
        dims_o(q, :) = [rg_yo(qy,2)-rg_yo(qy,1)+1, rg_xo(qx,2)-rg_xo(qx,1)+1];
    end
end

% total number of workers (Q: facets workers, K: data workers)
numworkers = Q + K;
cirrus_cluster = parcluster('local');
cirrus_cluster.NumWorkers = numworkers;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
if cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads > ncores
    exit(1);
end

parpool(cirrus_cluster, numworkers);

% define parallel constants (known by each worker)
Qyp = parallel.pool.Constant(Qy);
Qxp = parallel.pool.Constant(Qx);
Qp = parallel.pool.Constant(Q);
Kp = parallel.pool.Constant(K);
c_chunksp = parallel.pool.Constant(c_chunks);
waveletp = parallel.pool.Constant(wavelet);
nlevelp = parallel.pool.Constant(nlevel);
offsetp = parallel.pool.Constant(offset);

% define composite variables (local to a given worker)
% wavelet auxiliary variables
% /!\ only simple indexing allowed into Composite objects from the master
% node
Iq = Composite();
dims_q = Composite();
temLIdxs_q = Composite();
temRIdxs_q = Composite();
I_overlap_q = Composite();
dims_overlap_q = Composite();
dims_overlap_ref_q = Composite();
status_q = Composite();
Ncoefs_q = Composite();
offsetLq = Composite();
offsetRq = Composite();
% dimension of the ghost cells
overlap_g_south = Composite();
overlap_g_east = Composite();
overlap_g_south_east = Composite();
overlap = Composite();
% constant overlap: facet size
dims_oq = Composite();
crop_nuclear = Composite();
crop_l21 = Composite();

% initialize composite variables and constants
for q = 1:Q
    Iq{q} = I(q, :);
    dims_q{q} = dims(q, :);
    temLIdxs_q{q} = temLIdxs{q};
    temRIdxs_q{q} = temRIdxs{q};
    I_overlap_q{q} = I_overlap{q};
    dims_overlap_q{q} = dims_overlap{q};
    status_q{q} = status(q, :);
    Ncoefs_q{q} = Ncoefs{q};
    
    % additional composite variables (for the zero padding, see if fewer elements can be used)
    dims_overlap_ref_q{q} = dims_overlap_ref(q,:);
    offsetLq{q} = offsetL(q,:);
    offsetRq{q} = offsetR(q,:);
    
    % check which overlap is the largest (sdwt2 or nuclear norms)
    bool_crop = dims_overlap_ref(q,:) >= dims_o(q,:); % 0: nuclear norm largest, 1: dwt2 largest
    tmp_crop_l21 = [0, 0];
    tmp_crop_nuclear = [0, 0];
    if bool_crop(1)
        % sdwt2 largest
        tmp_crop_nuclear(1) = dims_overlap_ref(q,1) - dims_o(q,1);
    else
        tmp_crop_l21(1) = dims_o(q,1) - dims_overlap_ref(q,1);
    end
    
    if bool_crop(2)
        % sdwt2 largest
        tmp_crop_nuclear(2) = dims_overlap_ref(q,2) - dims_o(q,2);
    else
        tmp_crop_l21(2) = dims_o(q,2) - dims_overlap_ref(q,2);
    end
    crop_l21{q} = tmp_crop_l21;
    crop_nuclear{q} = tmp_crop_nuclear;    
    overlap{q} = max(max(dims_overlap{q}) - dims(q,:), dims_o(q, :) - dims(q,:)); 
    % amount of overlap necessary for each facet
    dims_oq{q} = dims_o(q, :);
end

% overlap dimension of the neighbour (necessary to define the ghost cells properly)
% define weights, depending on the weigthing option
if strcmp(window_type, 'pc') % 'piecewise_constant'
    % create weight matrix Wo (if needed)
    Wo = zeros(M, N);
    for q = 1:Q
       Wo(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2)) = ...
           Wo(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2)) + ones(dims_o(q,:)); 
    end
    Wo = 1./Wo;
end

w = Composite();
for q = 1:Q
    [qy, qx] = ind2sub([Qy, Qx], q);
    if qy < Qy
        % S (qy+1, qx)
        overlap_g_south{q} = overlap{(qx-1)*Qy + qy+1};
        
        if qx < Qx
            % SE (qy+1, qx+1)
            overlap_g_south_east{q} = overlap{qx*Qy + qy+1};
        else
            overlap_g_south_east{q} = [0, 0];
        end
    else
        overlap_g_south{q} = [0, 0];
        overlap_g_south_east{q} = [0, 0];
    end
    if qx < Qx
        % E (qy, qx+1)
        overlap_g_east{q} = overlap{qx*Qy + qy};
    else
        overlap_g_east{q} = [0, 0];
    end
    
    % define the weigths (depends on the position of the facet inside the grid)
    % to be tested on an independent script first
    switch window_type
        case 'triangular'
            tol = 1e-3;
            if qx == 1
                wdx = [ones(1, dims(q,2)-d), linspace(1-tol, tol, d)];
            elseif qx == Qx
                wdx = [linspace(tol, 1-tol, d), ones(1, dims_o(q,2)-2*d), ones(1, d)];
            else
                wdx = [linspace(tol, 1-tol, d), ones(1, dims_o(q,2)-2*d), linspace(1-tol, tol, d)];
            end    
            if qy == 1
                wdy = [ones(1, dims(q,1)-d), linspace(1-tol, tol, d)];
            elseif qy == Qy
                wdy = [linspace(tol, 1-tol, d), ones(1, dims_o(q,1)-2*d), ones(1, d)];
            else
                wdy = [linspace(tol, 1-tol, d), ones(1, dims_o(q,1)-2*d), linspace(1-tol, tol, d)];
            end
            w{q} = (wdy.').*wdx; 
        case 'hamming'
            dims_diff = dims_o(q, :) - dims(q, :);
            if qx > 1
                wc = window('hamming',2*dims_diff(2)).';
                wc = [wc(1:dims_diff(2)), ones(1, dims(q, 2))];
            else
                wc = ones(1, dims_o(q, 2));
            end
            if qy > 1
                wr = window('hamming',2*dims_diff(1)).';
                wr =  [wr(1:dims_diff(1)), ones(1, dims(q, 1))];
            else
                wr = ones(1, dims_o(q, 1));
            end
            w{q} = (wr.').*wc;
        case 'pc' % 'piecewise_constant'
            w{q} = Wo(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2));
        otherwise % make sure there is no 0 at the boundaries of the window
            w{q} = ones(dims_o(q,:));
    end 
end
clear XX YY xx yy Xq Yq V Wo
%%-- end initialisation auxiliary variables sdwt2

% Initializations.
init_flag = isfile(init_file_name);
if init_flag
    init_m = matfile(init_file_name);
end

if init_flag
    xsol = init_m.xsol;
    param = init_m.param;
    epsilon = init_m.epsilon;
    fprintf('xsol, param and epsilon uploaded \n\n')
else
    xsol = zeros(M,N,c);
    fprintf('xsol initialized \n\n')
end
% Primal / prior nodes (l21/nuclear norm dual variables)
v0_ = Composite();
weights0_ = Composite();
v1_ = Composite();
weights1_ = Composite();
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
        end
    end
    fprintf('v0, v1, weigths0, weights1 uploaded \n\n')
else
    spmd
        if labindex <= Qp.Value
            max_dims = max(dims_overlap_ref_q, dims_oq);
            [v0_, v1_, weights0_, weights1_] = initialize_dual_variables_prior_cst_overlap(Ncoefs_q, max_dims-crop_nuclear, c, nlevelp.Value);
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
adapt_eps_rel_var = parallel.pool.Constant(param.adapt_eps_rel_var);
adapt_eps_change_percentage = parallel.pool.Constant(param.adapt_eps_change_percentage);

Ap = Composite();
Atp = Composite();
Hp = Composite();
xhat_i = Composite();
Tp = Composite();
yp = Composite();
pUp = Composite();
Wp = Composite();
epsilonp = Composite();

% [09/10/2019] fixing bug in initialization of norm_res (warm-restart)
norm_res = Composite();
if init_flag
    for k = 1:K
        norm_res(Q+k) = init_m.norm_res(k,1);
    end
else
    for k = 1:K
        norm_res_tmp = cell(length(c_chunks{k}), 1);
        for i = 1:length(c_chunks{k})
            norm_res_tmp{i} = cell(length(y{k}{i}),1);
            for j = 1 : length(y{k}{i})
                norm_res_tmp{i}{j} = norm(y{k}{i}{j});
            end
        end
        norm_res{Q+k} = norm_res_tmp;
    end
    clear norm_res_tmp
end

sz_y = cell(K, 1);
n_blocks = cell(K, 1);
for k = 1:K
    sz_y{k} = cell(length(c_chunks{k}), 1);
    n_blocks{k} = cell(length(c_chunks{k}), 1);
    for i = 1:length(c_chunks{k})
        n_blocks{k}{i} = length(y{k}{i});
        for j = 1 : length(y{k}{i})
            sz_y{k}{i}{j} = numel(y{k}{i}{j});
        end
    end
    yp{Q+k} = y{k};
    y{k} = [];
    xhat_i{Q+k} = zeros(M, N, length(c_chunks{k}));
    Ap{Q+k} = A;
    Atp{Q+k} = At;
    Hp{Q+k} = H{k};
    Tp{Q+k} = T{k};
    T{k} = [];
    Wp{Q+k} = W{k};
    pUp{Q+k} = pU{k};
    epsilonp{Q+k} = epsilon{k};
end
clear epsilon pU W A At H y

v2_ = Composite();
t_block = Composite();
proj_ = Composite();
if init_flag
    for k = 1:K
        v2_(Q+k) = init_m.v2(k,1);
        proj_(Q+k) = init_m.proj(k,1);
        t_block(Q+k) = init_m.t_block(k,1);
    end
    fprintf('v2, proj, t_block uploaded \n\n')
else
    for k = 1:K
        v2_tmp = cell(length(c_chunks{k}), 1);
        t_block_ = cell(length(c_chunks{k}), 1);
        proj_tmp = cell(length(c_chunks{k}), 1);
        for i = 1:length(c_chunks{k})
            v2_tmp{i} = cell(n_blocks{k}{i},1);
            t_block_{i} = cell(n_blocks{k}{i},1);
            proj_tmp{i} = cell(n_blocks{k}{i},1);
            for j = 1 : n_blocks{k}{i}
                v2_tmp{i}{j} = zeros(sz_y{k}{i}{j} ,1);
                t_block_{i}{j} = 0;
                proj_tmp{i}{j} = zeros(sz_y{k}{i}{j}, 1);
            end
        end
        v2_{Q+k} = v2_tmp;
        proj_{Q+k} = proj_tmp;
        t_block{Q+k} = t_block_;
    end
    fprintf('v2, proj, t_block initialized \n\n')
    clear proj_tmp v2_tmp t_block_
end

if isfield(param,'init_reweight_step_count')
    reweight_step_count = param.init_reweight_step_count;
    fprintf('reweight_step_count uploaded')
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
rw_counts = 1;


%% Reweighting parameters
sig_bar = param.sig_bar;
sig = param.sig;
reweight_alpha = param.reweight_alpha;
reweight_alphap = Composite();
for q = 1:Q
    reweight_alphap{q} = reweight_alpha;
end
reweight_alpha_ffp = parallel.pool.Constant(param.reweight_alpha_ff);
reweight_steps = param.reweight_steps;

g_q = Composite();
xsol_q = Composite();
if init_flag
    for q = 1:Q
        xsol_q{q} = xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :);
        g_q{q} = init_m.g(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :);
    end
    fprintf('g uploaded \n\n')
else
    for q = 1:Q
        xsol_q{q} = xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :);
        g_q{q} = zeros([dims(q, :), c]);
    end
    fprintf('g initialized \n\n')
end

%Step sizes computation
%Step size for the dual variables
sigma0 = 1.0/param.nu0;
sigma1 = 1.0/param.nu1;
sigma2 = 1.0/param.nu2;

% Step size primal
tau = 0.99/(sigma0*param.nu0 + sigma1*param.nu1 + sigma2*param.nu2);

% Update constant dual variables
sigma00 = parallel.pool.Constant(tau*sigma0);
sigma11 = parallel.pool.Constant(tau*sigma1);
sigma22 = parallel.pool.Constant(tau*sigma2);
beta0 = parallel.pool.Constant(param.gamma0/sigma0); % only needed on the "primal/prior" workers
beta1 = parallel.pool.Constant(param.gamma/sigma1);

% Variables for the stopping criterion
flag = 0;
rel_val = zeros(param.max_iter, 1);
end_iter = zeros(param.max_iter, 1);

if isfield(param, 'init_t_start')
    t_start = param.init_t_start;
    fprintf('t_start uploaded \n\n')
else
    param.init_t_start = 1;
    t_start = 1;
    fprintf('t_start initialized \n\n')
end

start_loop = tic;

for t = t_start : param.max_iter
    
    %fprintf('Iter %i\n',t);
    start_iter = tic;
    
    spmd
        if labindex <= Q
            % primal/prior nodes (1:Q)
            
            % update primal variable
            [xsol_q, xhat_q, rel_x_q, norm_x_q] = update_primal(xsol_q, g_q);
            
            % send xhat_q (communication towards the data nodes)
            for i = 1:K
                labSend(xhat_q(:,:,c_chunksp.Value{i}), Qp.Value+i);
            end
            
            % update ghost cells (versions of xhat with overlap)
            x_overlap = zeros([max_dims, size(xsol_q, 3)]);
            x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xhat_q;
            x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
            
            % update dual variables (nuclear, l21) % errors here
            [v0_, g0] = update_nuclear_spmd_weighted(v0_, x_overlap(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :), w, weights0_, beta0.Value);
            [v1_, g1] = update_l21_spmd(v1_, x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :), weights1_, beta1.Value, Iq, ...
                dims_q, I_overlap_q, dims_overlap_q, offsetp.Value, status_q, ...
                nlevelp.Value, waveletp.Value, Ncoefs_q, temLIdxs_q, temRIdxs_q, offsetLq, offsetRq, dims_overlap_ref_q);
            g = zeros(size(x_overlap));
            g(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :) = sigma00.Value*g0;
            g(crop_l21(1)+1:end, crop_l21(2)+1:end, :) = g(crop_l21(1)+1:end, crop_l21(2)+1:end, :) + sigma11.Value*g1;          
            g = comm2d_reduce(g, overlap, Qyp.Value, Qxp.Value);
            
            % compute g_ for the final update term
            g_q = g(overlap(1)+1:end, overlap(2)+1:end, :);
            
            % retrieve portions of g2 from the data nodes
            for i = 1:Kp.Value
                g_q(:,:,c_chunksp.Value{i}) = g_q(:,:,c_chunksp.Value{i}) + labReceive(Qp.Value+i);
            end
        else
            % data nodes (Q+1:Q+K) (no data blocking, just frequency for
            % the moment)
            % retrieve xhat_i from the prior/primal nodes
            for q = 1:Qp.Value
                xhat_i(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2), :) = ...
                    labReceive(q);
            end
            [v2_, g2, proj_, norm_res, norm_residual_check_i, norm_epsilon_check_i] = update_data_fidelity_dr_block(v2_, yp, xhat_i, proj_, Ap, Atp, Hp, Tp, Wp, pUp, epsilonp, ...
                elipse_proj_max_iter.Value, elipse_proj_min_iter.Value, elipse_proj_eps.Value, sigma22.Value); % *_dr version when no blocking
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
    fprintf('Iter = %i, Time = %e\n',t,end_iter(t));
    
    %% Display
    if ~mod(t,500)
        
        %% compute value of the priors in parallel
        spmd
            if labindex <= Qp.Value
                x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
                x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
                [l21_norm, nuclear_norm] = prior_overlap_spmd_cst3_weighted(x_overlap, Iq, ...
                    offset, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
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
        
        % retrieve value of the monitoring variables (residual norms + epsilons)
        norm_epsilon_check = 0;
        norm_residual_check = 0;
        for i = Q+1:Q+K
            norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{i};
            norm_residual_check = norm_residual_check + norm_residual_check_i{i};
        end
        norm_epsilon_check = sqrt(norm_epsilon_check);
        norm_residual_check = sqrt(norm_residual_check);
        %% --
        
        % Log
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
            fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
        end
        
%         fitswrite(x0, ['./x0_it', num2str(t), '.fits']);
%         fitswrite(xsol, ['./xsol_it', num2str(t), '.fits']);
    end
    
    %% Global stopping criteria
    % if t>1 && rel_val(t) < param.rel_var && reweight_step_count > param.total_reweights && ...
    %         (norm_residual_check <= param.adapt_eps_tol_out*norm_epsilon_check)
    if ((t>1) && (reweight_step_count >= param.total_reweights)) && ((rel_val(t) < param.rel_var && ...
        (norm(residual_check) < param.adapt_eps_tol_out*norm(epsilon_check))) || ...
        (t - reweight_last_step_iter >= param.ppd_max_iter))
        flag = 1;
        break;
    end
    
    %% Update epsilons (in parallel)
    if param.use_adapt_eps && (t > param.adapt_eps_start) && (rel_val(t) < param.adapt_eps_rel_var)
        spmd
            if labindex > Qp.Value
                [epsilonp, t_block] = update_epsilon(epsilonp, t, t_block, rel_val(t), norm_res, ...
                    adapt_eps_tol_in.Value, adapt_eps_tol_out.Value, adapt_eps_steps.Value, adapt_eps_rel_var.Value, ...
                    adapt_eps_change_percentage.Value);
            end
        end
    end
    
    %% Reweighting (in parallel)
    % if (param.use_reweight_steps && (rel_val(t) < param.reweight_rel_var) && ...
    %     (reweight_step_count <= param.total_reweights) && ...
    %     (norm_residual_check <= param.adapt_eps_tol_out*norm_epsilon_check))
    is_converged_ppd = ((t - reweight_last_step_iter) >= param.ppd_min_iter) && (((rel_val(t) <= param.reweight_rel_var) && ...
(norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check))) || ...
((t - reweight_last_step_iter) >= param.ppd_max_iter));
        
    if is_converged_ppd && (reweight_step_count < param.total_reweights) % corresponds to the PPD stopping criterion
        fprintf('Reweighting: %i\n\n', reweight_step_count);
        
        % get xsol back from the workers
        for q = 1:Q
            xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
        end
        
        spmd
            if labindex <= Qp.Value
                % update weights
                x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
                x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);

                [weights1_, weights0_] = update_weights_overlap2_weighted(x_overlap, size(v1_), ...
                    Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, ...
                    Ncoefs_q, dims_overlap_ref_q, offsetLq, offsetRq, ...
                    reweight_alphap, crop_l21, crop_nuclear, w, sig, sig_bar);
                reweight_alphap = max(reweight_alpha_ffp.Value*reweight_alphap, 1);
            else
                % compute residual image on the data nodes
                res_ = compute_residual_images_dr_block(xsol(:,:,c_chunks{labindex-Qp.Value}), yp, Tp, Ap, Atp, Hp, Wp); % *_dr w/o data blocking
            end
        end
        reweight_alpha = max(param.reweight_alpha_ff*reweight_alpha, 1);     
        param.reweight_alpha = reweight_alpha;
        param.init_reweight_step_count = reweight_step_count+1;
        param.init_reweight_last_iter_step = t;
        param.init_t_start = t;
        
        if (reweight_step_count == 0) || (reweight_step_count == 1) || (~mod(reweight_step_count,5))
            % Save parameters (matfile solution)
            m = matfile([name, '_', ...
              num2str(param.ind) '_' num2str(param.gamma) '_' num2str(reweight_step_count) '.mat'], ...
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
                m.res(:,:,c_chunks{k}) = res_{Q+k};
                res_{Q+k} = [];
                m.v2(k,1) = v2_(Q+k);
                m.proj(k,1) = proj_(Q+k);
                m.t_block(k,1) = t_block(Q+k);
                m.epsilon(k,1) = epsilonp(Q+k);
                m.norm_res(k,1) = norm_res(Q+k);
            end
            clear m
        end 
        
        if (reweight_step_count >= param.total_reweights)
            param.reweight_max_reweight_itr = t+1;
            fprintf('\n\n No more reweights \n\n');
            break;
        end
        
        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;
        rw_counts = rw_counts + 1;
    end
end
toc(start_loop)

% Collect image facets back to the master
for q = 1:Q
    xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
end

% [09/10/2019] Modification from previous version: return fewer parameters,
% save through matfile to reduce memory footprint.
% Calculate residual images
spmd
    if labindex <= Qp.Value
        x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
        x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
        
        [l21_norm, nuclear_norm] = prior_overlap_spmd_cst3_weighted(x_overlap, Iq, ...
            offset, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
            offsetLq, offsetRq, crop_l21, crop_nuclear, w, size(v1_));
    else
        res_ = compute_residual_images_dr_block(xsol(:,:,c_chunks{labindex-Qp.Value}), yp, Tp, Ap, Atp, Hp, Wp); % *_dr w/o data blocking
    end
end

m = matfile([name, '_', ...
              num2str(param.ind) '_' num2str(param.gamma) '_' num2str(reweight_step_count) '.mat'], ...
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
    m.res(:,:,c_chunks{k}) = res_{Q+k};
    res_{Q+k} = [];
    m.v2(k,1) = v2_(Q+k);
    v2_{Q+k} = [];
    m.proj(k,1) = proj_(Q+k);
    proj_{Q+k} = [];
    m.t_block(k,1) = t_block(Q+k);
    t_block{Q+k} = [];
    m.epsilon(k,1) = epsilonp(Q+k);
    epsilonp{Q+k} = [];
    m.norm_res(k,1) = norm_res(Q+k);
end
m.xsol = xsol;
epsilon = m.epsilon; % see if necessary
norm_res_out = sqrt(sum(sum(sum((m.res).^2))));

% Update param structure and save
param.reweight_alpha = reweight_alpha;
param.init_reweight_step_count = reweight_step_count;
param.init_reweight_last_iter_step = t;
param.init_t_start = t;
m.param = param;
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
for i = Q+1:Q+K
    norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{i};
    norm_residual_check = norm_residual_check + norm_residual_check_i{i};
end
norm_epsilon_check = sqrt(norm_epsilon_check);
norm_residual_check = sqrt(norm_residual_check);

if (param.verbose > 0)
    if (flag == 1)
        fprintf('Solution found\n');
        fprintf('Iter %i\n',t);
        fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
        fprintf('epsilon = %e, residual = %e\n', norm_epsilon_check,norm_residual_check);
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf('Iter %i\n',t);
        fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
        fprintf('epsilon = %e, residual = %e\n', norm_epsilon_check,norm_residual_check);
    end
end

% end_iter = end_iter(end_iter > 0);
% rel_val = rel_val(1:numel(end_iter));

end
