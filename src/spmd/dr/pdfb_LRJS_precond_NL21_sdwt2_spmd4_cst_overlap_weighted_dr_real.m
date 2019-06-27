function [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,res,end_iter] = ...
    pdfb_LRJS_precond_NL21_sdwt2_spmd4_cst_overlap_weighted_dr_real(y, imsize, ...
    epsilon, A, At, H, pU, T, W, param, Qx, Qy, K, wavelet, L, nlevel, c_chunks, c, d, window_type)

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
if strcmp(window_type, 'piecewise_constant')
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
        case 'piecewise_constant'
            w{q} = Wo(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2));
        otherwise % make sure there is no 0 at the boundaries of the window
            w{q} = ones(dims_o(q,:));
    end 
end
clear XX YY xx yy Xq Yq V Wo
%%-- end initialisation auxiliary variables sdwt2

%Initializations.
if isfield(param,'init_xsol')
    xsol = param.init_xsol;
    fprintf('xsol uploaded \n\n')
else
    xsol = zeros(M,N,c);
    param.init_xsol = xsol;
    fprintf('xsol initialized \n\n')
end

% Primal / prior nodes (l21/nuclear norm dual variables)
v0_ = Composite();
weights0_ = Composite();
v1_ = Composite();
weights1_ = Composite();
if isfield(param,'init_v0') || isfield(param,'init_v1')
    for q = 1:Q
        v0_{q} = param.init_v0{q};
        v1_{q} = param.init_v1{q};
        weights0_{q} = param.init_weights0{q};
        weights1_{q} = param.init_weights1{q};
        spmd
            if labindex <= Qp.Value
                max_dims = max(dims_overlap_ref_q, dims_oq);
            end
        end
    end
    fprintf('v0, v1, weigths0, weights1 uploaded \n\n')
else
    param.init_v0 = cell(Q, 1);
    param.init_v1 = cell(Q, 1);
    param.init_weights0 = cell(Q, 1);
    param.init_weights1 = cell(Q, 1);
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
adapt_eps_rel_obj = parallel.pool.Constant(param.adapt_eps_rel_obj);
adapt_eps_change_percentage = parallel.pool.Constant(param.adapt_eps_change_percentage);

Ap = Composite();
Atp = Composite();
Hp = Composite();
x_hat_i = Composite();
Tp = Composite();
yp = Composite();
pUp = Composite();
Wp = Composite();
epsilonp = Composite();
norm_res = Composite();
sz_y = cell(K, 1);
n_blocks = cell(K, 1);
for k = 1:K
    norm_res_tmp = cell(length(c_chunks{k}), 1);
    sz_y{k} = cell(length(c_chunks{k}), 1);
    n_blocks{k} = cell(length(c_chunks{k}), 1);
    for i = 1:length(c_chunks{k})
        norm_res_tmp{i} = cell(length(y{k}{i}),1);
        n_blocks{k}{i} = length(y{k}{i});
        for j = 1 : length(y{k}{i})
            norm_res_tmp{i}{j} = norm(y{k}{i}{j});
            sz_y{k}{i}{j} = numel(y{k}{i}{j});
        end
    end
    yp{Q+k} = y{k};
    y{k} = [];
    x_hat_i{Q+k} = zeros(M, N, length(c_chunks{k}));
    Ap{Q+k} = A;
    Atp{Q+k} = At;
    Hp{Q+k} = H(c_chunks{k});
    Tp{Q+k} = T{k};
    T{k} = [];
    Wp{Q+k} = W{k};
    pUp{Q+k} = pU{k};
    epsilonp{Q+k} = epsilon{k};
    norm_res{Q+k} = norm_res_tmp;
end

clear norm_res_tmp epsilon pU W A At H y

v2_ = Composite();
t_block = Composite();
proj_ = Composite();
if isfield(param,'init_v2') % assume all the other related elements are also available in this case
    for k = 1:K
        v2_{Q+k} = param.init_v2{k};
        proj_{Q+k} = param.init_proj{k};
        t_block{Q+k} = param.init_t_block{k};
    end
    fprintf('v2, proj, t_block uploaded \n\n')
else
    param.init_v2 = cell(K, 1);
    param.init_proj = cell(K, 1);
    param.init_t_block = cell(K, 1);
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
end

clear proj_tmp v2_tmp norm_res_tmp t_block_ G y

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

reweight_alpha = param.reweight_alpha;
reweight_alphap = Composite();
for q = 1:Q
    reweight_alphap{q} = reweight_alpha;
end
reweight_alpha_ffp = parallel.pool.Constant(param.reweight_alpha_ff);
reweight_steps = param.reweight_steps;

g_q = Composite();
xsol_q = Composite();
if isfield(param,'init_g')
    for q = 1:Q
        xsol_q{q} = xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :);
        g_q{q} = param.init_g(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :);
    end
    fprintf('g uploaded \n\n')
else
    param.init_g = zeros(size(xsol)); % 0 values...
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
rel_fval = zeros(param.max_iter, 1);
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
            [v2_, g2, proj_, norm_res, norm_residual_check_i, norm_epsilon_check_i] = update_data_fidelity_dr(v2_, yp, xhat_i, proj_, Ap, Atp, Hp, Tp, Wp, pUp, epsilonp, ...
                elipse_proj_max_iter.Value, elipse_proj_min_iter.Value, elipse_proj_eps.Value, sigma22.Value);
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
    rel_fval(t) = sqrt(rel_x/norm_x);
    end_iter(t) = toc(start_iter);
    fprintf('Iter = %i, Time = %e\n',t,end_iter(t));
    
    %% Display
    if ~mod(t,100)
        
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
            fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
            fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
        end
        
%         fitswrite(x0, './x0.fits');
%         fitswrite(xsol, './xsol.fits');
    end
    
    %% Global stopping criteria
    if t>1 && rel_fval(t) < param.rel_obj && reweight_step_count > param.total_reweights && ...
            (norm_residual_check <= param.adapt_eps_tol_out*norm_epsilon_check)
        flag = 1;
        break;
    end
    
    %% Update epsilons (in parallel)
    if param.use_adapt_eps && t > param.adapt_eps_start
        spmd
            if labindex > Qp.Value
                [epsilonp, t_block] = update_epsilon(epsilonp, t, t_block, rel_fval(t), norm_res, ...
                    adapt_eps_tol_in.Value, adapt_eps_tol_out.Value, adapt_eps_steps.Value, adapt_eps_rel_obj.Value, ...
                    adapt_eps_change_percentage.Value);
            end
        end
    end
    
    %% Reweighting (in parallel)
    if (param.step_flag && rel_fval(t) < param.reweight_rel_obj && ...
        norm_residual_check <= param.adapt_eps_tol_out*norm_epsilon_check)
        reweight_steps = (t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size));
        param.step_flag = 0;
    end
    
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
            (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
            t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);
        
        % get xsol back from the workers
        for q = 1:Q
            xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
        end
        
        % update weights
        spmd
            if labindex <= Qp.Value
                x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
                x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);

                [weights1_, weights0_] = update_weights_overlap2_weighted(x_overlap, size(v1_), ...
                    Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, ...
                    Ncoefs_q, dims_overlap_ref_q, offsetLq, offsetRq, ...
                    reweight_alphap, crop_l21, crop_nuclear, w);
                reweight_alphap = reweight_alpha_ffp.Value * reweight_alphap;
            else
                % compute residual image on the data nodes
                res_ = compute_residual_images_dr(xsol(:,:,c_chunks{labindex-Qp.Value}), yp, Tp, Ap, Atp, Hp, Wp);
            end
        end
        reweight_alpha = param.reweight_alpha_ff .* reweight_alpha; % on the master node       
        
        % Save parameters
        epsilon = cell(K, 1);
        for k = 1:K
            param.init_v2{k} = v2_{Q+k};
            param.init_proj{k} = proj_{Q+k};
            param.init_t_block{k} = t_block{Q+k};
            epsilon{k} = epsilonp{Q+k};
        end
        for q = 1:Q
            param.init_v0{q} = v0_{q};
            param.init_v1{q} = v1_{q};
            param.init_weights0{q} = weights0_{q};
            param.init_weights1{q} = weights1_{q};
            param.init_xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
            param.init_g(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = g_q{q};
        end
        param.reweight_alpha = reweight_alpha;
        param.init_reweight_step_count = reweight_step_count+1;
        param.init_reweight_last_iter_step = t;
        param.init_t_start = t;
        
        % Compute residual images (master node)
        res = zeros(size(xsol));
        for k = 1 : K
            res(:,:,c_chunks{k}) = res_{Q+k};
            res_{Q+k} = [];
        end
        
        if (reweight_step_count == 0) || (reweight_step_count == 1) || (~mod(reweight_step_count,5))
            mkdir('./results/')
            save(['./results/result_HyperSARA_spmd4_cst_weighted_rd_' num2str(param.ind) '_' num2str(param.gamma) '_' num2str(reweight_step_count) '.mat'],'-v7.3','param','res', 'epsilon');
        end 
        
        if (reweight_step_count >= param.total_reweights)
            param.reweight_max_reweight_itr = t+1;
            fprintf('\n\n No more reweights \n\n');
            break;
        end
        
        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;
        rw_counts = rw_counts + 1;
        
        epsilon = [];
    end
end
toc(start_loop)

% Collect image back to the master
for q = 1:Q
    xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
end

% Calculate residual images
res = zeros(size(xsol));
spmd
    if labindex > Qp.Value
        res_ = compute_residual_images_dr(xsol(:,:,c_chunks{labindex-Qp.Value}), yp, Tp, Ap, Atp, Hp, Wp);
    end
end
for k = 1 : K
    res(:,:,c_chunks{k}) = res_{Q+k};
end
norm_res_out = norm(res(:));

%Final log (merge this step with the computation of the residual image for
% each frequency of interest)
spmd
    if labindex <= Qp.Value
        x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
        x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
        
        [l21_norm, nuclear_norm] = prior_overlap_spmd_cst3_weighted(x_overlap, Iq, ...
            offset, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
            offsetLq, offsetRq, crop_l21, crop_nuclear, w, size(v1_));
    end
end

epsilon = cell(K, 1);
for k = 1:K
    param.init_v2{k} = v2_{Q+k};
    v2_{Q+k} = [];
    param.init_proj{k} = proj_{Q+k};
    proj_{Q+k} = [];
    param.init_t_block{k} = t_block{Q+k};
    t_block{Q+k} = [];
    epsilon{k} = epsilonp{Q+k};
    epsilonp{Q+k} = [];
end
for q = 1:Q
    param.init_v0{q} = v0_{q};
    v0_{q} = [];
    param.init_v1{q} = v1_{q}; 
    v1_{q} = [];
    param.init_weights0{q} = weights0_{q};
    weights0_{q} = [];
    param.init_weights1{q} = weights1_{q};
    weights1_{q} = [];
    param.init_xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
    xsol_q{q} = [];
    param.init_g(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = g_q{q};
    g_q{q} = [];
end
param.reweight_alpha = reweight_alpha;
param.init_reweight_step_count = reweight_step_count;
param.init_reweight_last_iter_step = t;
param.init_t_start = t;

l21 = 0;
nuclear = 0;
for q = 1:Q
    l21 = l21 + l21_norm{q};
    nuclear = nuclear + nuclear_norm{q};
end

if (param.verbose > 0)
    if (flag == 1)
        fprintf('Solution found\n');
        fprintf('Iter %i\n',t);
        fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
        fprintf('epsilon = %e, residual = %e\n', norm_epsilon_check,norm_residual_check);
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf('Iter %i\n',t);
        fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
        fprintf('epsilon = %e, residual = %e\n', norm_epsilon_check,norm_residual_check);
    end
end

% end_iter = end_iter(end_iter > 0);
% rel_fval = rel_fval(1:numel(end_iter));

end
