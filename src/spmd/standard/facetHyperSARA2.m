function [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,SNR,SNR_average] = ...
    facetHyperSARA2(y, epsilon, A, At, pU, G, W, param, X0, Qx, Qy, K, wavelet, ...
    L, nlevel, c_chunks, c, init_file_name)
%facetHyperSARA2: ...
%
% version with the same tessellation for both the facet nuclear and l21
% norms, no spatial weighting correction for the nuclear norm. Alternative implementation of
% facetHyperSARA: the image is updated and stored on the data nodes, 
% instead of the facet nodes.
%
%-------------------------------------------------------------------------%
%%
% Input: 
%
% > y           blocks of visivilities {L}{nblocks_l}
% > epsilon     l2-ball norms {L}{nblocks_l}
% > A           measurement operator
% > At          adjoint of the measurement operator
% > pU          preconditioning matrices {L}{nblocks_l}
% > G           gridding matrices {L}{nblocks_l}
% > W           masks for selection of the blocks of visibilities
% > param       algorithm parameters (struct)
%
%   general
%   > .verbose           print log or not
%   > .rel_obj   (1e-5)  stopping criterion
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
%   > .reweight_rel_obj   (1e-4)     criterion for performing reweighting
%   > .reweight_min_steps_rel_obj (300) minimum number of iterations between consecutive reweights
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
%   > .adapt_eps_rel_obj (5e-4)      bound on the relative change of the solution
%   > .adapt_eps_change_percentage  (0.5*(sqrt(5)-1)) weight of the update w.r.t the l2 norm of the residual data
%
%
% > X0          ground truth wideband image [L, N]
% > Qx          number of facets along dimension x [1]
% > Qy          number of facets along dimension y [1]
% > K           number of datab computing processes [1]
% > wavelet     wavelet doctionaries considered (should contain 'self' by
%               default in last position)
% > L           size of the wavelet filters considered (by cinvention, 0 for the Dirac basis)
% > nlevel      decomposition depth [1]
% > c_chunks    indices of the bands handled by each data node {K, 1}
% > c           total number of specrtal channels [1]
% > init_file_name  name of a valid .mat file for initialization (for warm-restart)
%
%
% Output:
%
% < xsol        reconstructed wideband image [L, N]
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
% < rel_fval        relative variation
% < nuclear         value of the faceted nuclear norm
% < l21             value of the l21 regularization term
% < norm_res_out    norm of the reidual image 
% < res             residual image [N(1), N(2)]
% < end_iter        last iteration
%-------------------------------------------------------------------------%
%%
% Data update and primal update gathered in the same nodes 

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
%%
%% TO BE FINALIZED, BUT DOES NOT MAKE SENSE FOR THE REWEIGTHING + COMPUTATION OF THE PRIORS IN PARALLEL...
% (EXTRA COMMUNICATIONS NEEDED)

% maxNumCompThreads(param.num_workers);

% initialize monitoring variables (display active)
SNR = 0;
SNR_average = 0;
norm_epsilon_check = 0;
norm_residual_check = 0;

% size of the oversampled Fourier space (vectorized)
No = size(W{1}{1}{1}, 1);

% number of pixels (spatial dimensions)
[M, N] = size(At(zeros(No, 1)));

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
[~, ~, ~, dims_overlap_ref, I_overlap, dims_overlap, ...
    ~, ~, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices([M, N], I, dims, nlevel, wavelet, L);
id_dirac = find(ismember(wavelet, 'self'), 1);
dirac_present = ~isempty(id_dirac);

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

% define parallel constants(known by each worker)
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
xhat_q = Composite();
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
    overlap{q} = max(dims_overlap{q}) - dims(q,:); % amount of overlap necessary for each facet
    xhat_q{q} = zeros([dims(q, :), c]);
end

% amount of overlap of the neighbour (necessary to define the ghost cells properly)
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
end
%%-- end initialisation auxiliary variables sdwt2

%Initializations.
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

% Prior/primal nodes
% l21 / nuclear norm dual variables
v0_ = Composite();
weights0_ = Composite();
v1_ = Composite();
weights1_ = Composite(); % idem, can be created in parallel with spmd
if init_flag
    for q = 1:Q
        v0_(q) = init_m.v0(q,1);
        v1_(q) = init_m.v1(q,1);
        weights0_(q) = init_m.weights0(q,1);
        weights1_(q) = init_m.weights1(q,1);
    end
else
    spmd
        if labindex <= Qp.Value
            [v0_, v1_, weights0_, weights1_] = initialize_dual_variables_prior_overlap(Ncoefs_q, dims_q, dims_overlap_ref_q, dirac_present, c, nlevelp.Value);
        end
    end
end

%% Data node parameters

A = afclean(A);
At = afclean(At);

% see if these constants are properly updated in practice
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
Gp = Composite();
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
    Ap{Q+k} = A;
    Atp{Q+k} = At;
    Gp{Q+k} = G{k};
    G{k} = [];
    Wp{Q+k} = W{k};
    pUp{Q+k} = pU{k};
    epsilonp{Q+k} = epsilon{k};
end
clear epsilon pU W G y

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
    clear proj_tmp v2_tmp t_block_
    fprintf('v2, proj, t_block initialized \n\n')
end

if isfield(param,'init_reweight_step_count')
    reweight_step_count = param.init_reweight_step_count;
    fprintf('reweight_step_count uploaded')
else
    reweight_step_count = 0;
    fprintf('reweight_step_count initialized \n\n')
end

if isfield(param,'init_reweight_last_iter_step')
    reweight_last_step_iter = param.init_reweight_last_iter_step;
    fprintf('reweight_last_iter_step uploaded \n\n')
else
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

g_i = Composite();
xsol_i = Composite();
if init_flag
    for i = 1:K
        xsol_i{Q+k} = xsol(:, :, c_chunks{k});
        g_i{Q+i} = init_m.g(:,:, c_chunks{k});
    end
    fprintf('g uploaded \n\n')
else
    for i = 1:K
        xsol_i{Q+k} = xsol(:, :, c_chunks{k});
        g_i{Q+i} = zeros([M, N, numel(c_chunks{i})]);
    end
    fprintf('g initialized \n\n')
end

%Step sizes computation
%Step size for the dual variables
sigma0 = 1.0/param.nu0;
sigma1 = 1.0/param.nu1;
sigma2 = 1.0/param.nu2;

%Step size primal
tau = 0.99/(sigma0*param.nu0 + sigma1*param.nu1 + sigma2*param.nu2);

% Update constant dual variables
sigma00 = parallel.pool.Constant(tau*sigma0);
sigma11 = parallel.pool.Constant(tau*sigma1);
sigma22 = parallel.pool.Constant(tau*sigma2);
beta0 = parallel.pool.Constant(param.gamma0/sigma0); % needed only on the "prior" workers
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
            % retrieve xhat_q from the primal nodes
            for i = 1:Kp.Value
                xhat_q(:, :, c_chunksp.Value{i}) = labReceive(Qp.Value+i);
            end
            
            % update ghost cells (versions of xhat with overlap)
            % overlap_q = dims_overlap_ref_q - dims_q;
            x_overlap = zeros([dims_overlap_ref_q, size(xhat_q, 3)]);
            x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xhat_q;
            x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
            
            % update dual variables (nuclear, l21)
            [v0_, g0] = update_nuclear_spmd(v0_, x_overlap, weights0_, beta0.Value);
            [v1_, g1] = update_l21_spmd(v1_, x_overlap, weights1_, beta1.Value, Iq, ...
                dims_q, I_overlap_q, dims_overlap_q, offsetp.Value, status_q, ...
                nlevelp.Value, waveletp.Value, Ncoefs_q, temLIdxs_q, temRIdxs_q, offsetLq, offsetRq, dims_overlap_ref_q);
            g = sigma00.Value*g0 + sigma11.Value*g1;
            g = comm2d_reduce(g, overlap, Qyp.Value, Qxp.Value);
            
            % compute g_ for the final update term
            g_q = g(overlap(1)+1:end, overlap(2)+1:end, :);
            
            % send portions of g_q to the data nodes
            for i = 1:Kp.Value
                labSend(g_q(:,:,c_chunksp.Value{i}), Qp.Value+i);
            end
        else
            % update primal variable
            [xsol_i, xhat_i, rel_x_i, norm_x_i] = update_primal(xsol_i, g_i);
            
            % send xhat_i (communication towards the data nodes)
            for q = 1:Qp.Value
                labSend(xsol_i(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2),:,:), q); % to be changed
            end
            
            % data nodes (Q+1:Q+K) (no data blocking, just frequency for the moment
            [v2_, g_i, proj_, norm_res, norm_residual_check_i, norm_epsilon_check_i] = update_data_fidelity(v2_, yp, xhat_i, proj_, Ap, Atp, Gp, Wp, pUp, epsilonp, ...
                elipse_proj_max_iter.Value, elipse_proj_min_iter.Value, elipse_proj_eps.Value, sigma22.Value);
            
            % retrieve portions of g2 to the prior/primal nodes
            for q = 1:Qp.Value
                g_i(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2), :) = ...
                g_i(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2), :) + labReceive(q); % to be done inside the function?
            end
        end
    end
    
    %% Relative change of objective function
    % retrieve rel_x_q, norm_x_q for the workers
    rel_x = 0;
    norm_x = 0;
    for i = 1:K
        rel_x = rel_x + rel_x_i{Q+i};
        norm_x = norm_x + norm_x_i{Q+i};
    end
    rel_fval(t) = sqrt(rel_x/norm_x);
    end_iter(t) = toc(start_iter);
    fprintf('Iter = %i, Time = %e\n',t,end_iter(t));
    
    %% Retrieve value of the monitoring variables (residual norms + epsilons)
    norm_epsilon_check = 0;
    norm_residual_check = 0;
    for i = Q+1:Q+K
        norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{i};
        norm_residual_check = norm_residual_check + norm_residual_check_i{i};
    end
    norm_epsilon_check = sqrt(norm_epsilon_check);
    norm_residual_check = sqrt(norm_residual_check);
    
    %% Display
    if ~mod(t,100)
        
        %% compute value of the priors in parallel
        spmd
            if labindex <= Qp.Value
                % compute values for the prior terms
                %x_overlap = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
                %see if this is necessary or not (to be further investigated)
                for i = 1:Kp.Value
                    xhat_q(:, :, c_chunksp.Value{i}) = labReceive(Qp.Value+i);
                end
            
                % update ghost cells (versions of xhat with overlap)
                % overlap_q = dims_overlap_ref_q - dims_q;
                x_overlap = zeros([dims_overlap_ref_q, size(xhat_q, 3)]);
                x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xhat_q;
                x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value); % problem index in l. 80 (position 2...) on worker 2 (to be investigated further)

                [l21_norm, nuclear_norm] = prior_overlap_spmd(x_overlap, Iq, ...
                    offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
                    offsetLq, offsetRq, size(v1_)); % ISSUE HERE...
            else
                for q = 1:Qp.Value
                    labSend(xsol_i(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2),:), q); % I, dims need to be known on each data node
                end
            end
        end
        
        % retrieve value of the priors
        l21 = 0;
        nuclear = 0;
        for q = 1:Q
            l21 = l21 + l21_norm{q};
            nuclear = nuclear + nuclear_norm{q};
        end
        
        % SNR
        % get xsol back from the workers
        for i = 1:K
            xsol(:, :, c_chunks{i}) = xsol_i{Q+i};
        end
        sol = reshape(xsol(:),numel(xsol(:))/c,c);
        SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
        psnrh = zeros(c,1);
        for i = 1:c
            psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
        end
        SNR_average = mean(psnrh);
        %% --
        
        % Log
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
            fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
            fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
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
    if (param.step_flag && t>500) % rel_fval(t) < param.reweight_rel_obj)
        reweight_steps = (t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size));
        param.step_flag = 0;
    end
    
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
            (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
            t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);

        spmd
            if labindex <= Qp.Value
                for i = 1:Kp.Value
                    xhat_q(:, :, c_chunksp.Value{i}) = labReceive(Qp.Value+i); % contains xsol here
                end
                
                x_overlap = zeros([dims_overlap_ref_q, size(xhat_q, 3)]);
                x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xhat_q;
                x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);

                [weights1_, weights0_] = update_weights_overlap(x_overlap, size(v1_), ...
                    Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, ...
                    Ncoefs_q, dims_overlap_ref_q, offsetLq, offsetRq, reweight_alphap);
                reweight_alphap = reweight_alpha_ffp.Value * reweight_alphap;
            else
                % compute residual image on the data nodes
                for q = 1:Qp.Value
                    labSend(xsol_i(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2),:), q); % I, dims need to be known on each data node
                end
                % compute residual image on the data nodes
                res_ = compute_residual_images(xsol_i, yp, Gp, Ap, Atp, Wp);
            end
        end
        reweight_alpha = param.reweight_alpha_ff .* reweight_alpha; % on the master node
        reweight_alpha = param.reweight_alpha_ff .* reweight_alpha; % on the master node
        param.reweight_alpha = reweight_alpha;
        param.init_reweight_step_count = reweight_step_count+1;
        param.init_reweight_last_iter_step = t;
        param.init_t_start = t; 
        %-- end modifications
        
        if (reweight_step_count == 0) || (reweight_step_count == 1) || (~mod(reweight_step_count,5))
            % Save parameters (matfile solution)
            mkdir('./results/')
            m = matfile(['./results/facetHyperSARA_' ...
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
            end
            % data nodes
            for k = 1:K
                m.xsol(:,:,c_chunks{k}) = xsol_i{Q+k};
                m.g(:,:,c_chunks{k}) = g_i{Q+k};
                m.res(:,:,c_chunks{k}) = res_{Q+k};
                res_{Q+k} = [];
                m.v2(k,1) = v2_(Q+k);
                m.proj(k,1) = proj_(Q+k);
                m.t_block(k,1) = t_block(Q+k);
                m.epsilon(k,1) = epsilonp(Q+k);
                m.norm_res(k,1) = norm_res(Q+k);
            end
            %m.SNR = SNR;
            %m.SNR_average = SNR_average;
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

spmd
    if labindex <= Qp.Value
        % compute values for the prior terms
        %x_overlap = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
        %see if this is necessary or not (to be further investigated)
        for i = 1:Kp.Value
            xhat_q(:, :, c_chunksp.Value{i}) = labReceive(Qp.Value+i); % contains xsol
        end
        
        % update ghost cells (versions of xhat with overlap)
        % overlap_q = dims_overlap_ref_q - dims_q;
        x_overlap = zeros([dims_overlap_ref_q, size(xhat_q, 3)]);
        x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xhat_q;
        x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
        
        [l21_norm, nuclear_norm] = prior_overlap_spmd(x_overlap, Iq, ...
            offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
            offsetLq, offsetRq, size(v1_));
    else
        for q = 1:Qp.Value
            labSend(xsol_i(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2),:), q);
        end
        % compute residual image on the data nodes
        res_ = compute_residual_images(xsol_i, yp, Gp, Ap, Atp, Wp);
    end
end

m = matfile(['./results/facetHyperSARA2_' ...
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
    m.xsol(:,:,c_chunks{k}) = xsol_i{Q+k};
    m.g(:,:,c_chunks{k}) = g_i{Q+k};
    g_i{Q+k} = [];
end
norm_res_out = sqrt(sum(sum(sum((m.res).^2))));
m.xsol = xsol;
epsilon = m.epsilon; % see if necessary

% Update param structure and save
param.reweight_alpha = reweight_alpha;
param.init_reweight_step_count = reweight_step_count;
param.init_reweight_last_iter_step = t;
param.init_t_start = t;
m.param = param;

% compute SNR (on the master node)
sol = reshape(xsol(:),numel(xsol(:))/c,c);
SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
psnrh = zeros(c,1);
for i = 1:c
    psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
end
SNR_average = mean(psnrh);
m.SNR = SNR;
m.SNR_average = SNR_average;
clear m

% % to be completely modified (within spmd function?)
% % Calculate residual images
% res = zeros(size(xsol));
% for k = 1 : K
%     for i = 1 : length(c_chunks{k})
%         Fx = A(xsol(:,:,c_chunks{k}(i)));
%         g2 = zeros(No,1);
%         for j = 1 : length(G{k}{i})
%             res_f = y{k}{i}{j} - G{k}{i}{j} * Fx(W{k}{i}{j});
%             u2 = G{k}{i}{j}' * res_f;
%             g2(W{k}{i}{j}) = g2(W{k}{i}{j}) + u2; % no overlap between different groups? (content of W...)
%         end
%         res(:,:,c_chunks{k}(i)) = real(At(g2)); % residual images
%     end
%     xsol(:, :, c_chunks{k}) = xsol_i{k};
% end

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
        fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
        fprintf('epsilon = %e, residual = %e\n', norm_epsilon_check,norm_residual_check);
        fprintf('SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf('Iter %i\n',t);
        fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
        fprintf('epsilon = %e, residual = %e\n', norm_epsilon_check,norm_residual_check);
        fprintf('SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
    end
end

end_iter = end_iter(end_iter > 0);
rel_fval = rel_fval(1:numel(end_iter));

end
