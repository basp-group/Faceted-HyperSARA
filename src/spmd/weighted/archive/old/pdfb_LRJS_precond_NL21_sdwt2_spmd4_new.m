function [xsol,v0,v1,v2,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = pdfb_LRJS_precond_NL21_sdwt2_spmd4_new(y, epsilon, A, At, pU, G, W, param, X0, Qx, Qy, K, wavelet, L, nlevel, c_chunks, c, p, bool_weights)

%SPMD version: use spmd for all the priors, deal with the data fidelity
% term in a single place.

% p: proportion of facet overlap
% 

%% NOTE:
% this version relies on a specialised version of sdwt2, slightly less
% general but faster (based on Arwa's work).

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
% 4. Backup options have been removed four the moment, see initialization
% from a set of known variables (warm-restart) [really useful in practice? 
%%

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

% define overlapping facets
rg_yo = domain_decomposition_overlap3(Qy, M, p);
rg_xo = domain_decomposition_overlap3(Qx, N, p);

% equivalent number of facets for the non-redundant version
Qx = floor(Qx/p)+1-floor(1/p);
Qy = floor(Qy/p)+1-floor(1/p);
Q = Qx*Qy;

% define starting indices/sizes of the overlapping facets 
Io = zeros(Q, 2);
dims_o = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        Io(q, :) = [rg_yo(qy, 1)-1, rg_xo(qx, 1)-1];
        dims_o(q, :) = [rg_yo(qy,2)-rg_yo(qy,1)+1, rg_xo(qx,2)-rg_xo(qx,1)+1];
    end
end

% create global weigth matrix (for l21 regularization term)

if bool_weights
    Weight = zeros(M, N);
    for q = 1:Q
       Weight(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2)) = ...
           Weight(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2)) + ones(dims_o(q,:)); 
    end
    Weight = 1./Weight;
else
    Weight = ones(M, N);
end

% define reference spatial facets (no overlap)
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
dims_oq = Composite();
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
% 
crop_nuclear = Composite();
crop_l21 = Composite();

% initialize composite variables and constants
for q = 1:Q
    Iq{q} = I(q,:);
    dims_q{q} = dims(q,:);
    dims_oq{q} = dims_o(q,:);
    temLIdxs_q{q} = temLIdxs{q};
    temRIdxs_q{q} = temRIdxs{q};
    I_overlap_q{q} = I_overlap{q};
    dims_overlap_q{q} = dims_overlap{q};
    status_q{q} = status(q,:);
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
%     if bool_crop{q}
%         overlap{q} = max(dims_overlap{q}) - dims(q,:); % amount of overlap necessary for each facet
%     else
%         overlap{q} = dims_o(q, :) - dims(q,:);
%     end
end

% overlap dimension of the neighbour (necessary to define the ghost cells properly)
%crop = Composite();
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
    
    % problem with crop: make sure I am in the apporpriate configuration!
    %crop{q} = abs(dims_o(q,:) - dims_overlap_ref(q,:)); % dims_overlap_ref(q,:) - dims_o(q,:);
    % assuming dims_o > dims_overlap_ref
    w{q} = Weight(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2));
end
clear qy qx;
%%-- end initialisation auxiliary variables sdwt2

%Initializations.
if isfield(param,'init_xsol')
    xsol = param.init_xsol;
    fprintf('xsol uploaded \n\n')
else
    xsol = zeros(M,N,c);
    fprintf('xsol NOT uploaded \n\n')
end

% Primal / prior nodes (l21/nuclear norm dual variables)
v0_ = Composite();
weights0_ = Composite();
v1_ = Composite();
weights1_ = Composite();
if isfield(param,'init_v0') || isfield(param,'init_v1')
    for q = 1:Q
        v0_{q} = param.init_v0{q};
        v1_{q} = param.init_v0{q};
        weights0_{q} = param.init_weights0{q};
        weights1_{q} = param.init_weights1{q};
    end
else
    spmd
        % to be changed (modify the sizes if needed, ans the way the elements are cropped)
        if labindex <= Qp.Value
            max_dims = max(dims_overlap_ref_q, dims_oq);
            [v0_, v1_, weights0_, weights1_] = initialize_dual_variables_prior_overlap(Ncoefs_q, dims_q, max_dims-crop_nuclear, dirac_present, c, nlevelp.Value);
            [qy, qx] = ind2sub([Qy, Qx], labindex); % error here... why (issue, only zeros for the image)
        end
    end
end

%% Data node parameters
A = afclean(A);
At = afclean(At);
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
x_hat_i = Composite();
Gp = Composite();
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
    Gp{Q+k} = G{k};
    G{k} = [];
    Wp{Q+k} = W{k};
    pUp{Q+k} = pU{k};
    epsilonp{Q+k} = epsilon{k};
    norm_res{Q+k} = norm_res_tmp;
end

clear norm_res_tmp epsilon pU W G y % W, G, y needed to compute the residual image 

v2_ = Composite();
t_block = Composite();
proj_ = Composite();
if isfield(param,'init_v2') % assume all the other related elements are also available in this case
    for k = 1:K
        v2_{Q+k} = param.init_v2{k};
        proj_{Q+k} = param.init_proj{k};
        t_block{Q+k} = param.init_t_block{k};
    end
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
end

clear proj_tmp v2_tmp norm_res_tmp t_block_ G y

reweight_last_step_iter = 0;
reweight_step_count = 0;
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
for q = 1:Q
    xsol_q{q} = xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :);
    g_q{q} = zeros([dims(q, :), c]);
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
t_start = 1; % use of t_start?

start_loop = tic;

for t = t_start : param.max_iter
    
    %fprintf('Iter %i\n',t);
    start_iter = tic;
    
    % issue at this stage (indices, to be investigated)
    spmd
        if labindex <= Q
            % primal/prior nodes (1:Q)
            
            % update primal variable
            [xsol_q, xhat_q, rel_x_q, norm_x_q] = update_primal(xsol_q, g_q);
            
            % force image to be zero on the borders (size of one overlap)
%             if (qx == 1) && (Qxp.Value > 1)
%                 xsol_q(:, 1:floor(p*dims_o(q,2)), :) = 0;
%                 xhat_q(:, 1:floor(p*dims_o(q,2)), :) = 0;
%             end
%             if (qx == Qxp.Value) && (Qxp.Value > 1)
%                 xsol_q(:, end-floor(p*dims_o(q,2))+1:end, :) = 0;
%                 xhat_q(:, end-floor(p*dims_o(q,2))+1:end, :) = 0;
%             end           
%             if (qy == 1) && (Qyp.Value > 1)
%                 xsol_q(1:floor(p*dims_o(q,1)), :, :) = 0;
%                 xhat_q(1:floor(p*dims_o(q,1)), :, :) = 0;
%             end
%             if (qy == Qyp.Value) && (Qyp.Value > 1)
%                 xsol_q(end-floor(p*dims_o(q,1))+1:end, :, :) = 0;
%                 xhat_q(end-floor(p*dims_o(q,1))+1:end, :, :) = 0;
%             end
      
            % send xhat_q (communication towards the data nodes)
            for i = 1:K
                labSend(xhat_q(:,:,c_chunksp.Value{i}), Qp.Value+i);
            end
            
            % update ghost cells (versions of xhat with overlap)
            x_overlap = zeros([max_dims, size(xsol_q, 3)]); % zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
            x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xhat_q;
            x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
            
            % update dual variables (nuclear, l21)
%             [v0_, g0] = update_nuclear_spmd_weighted(v0_, x_overlap(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :), w, weights0_, beta0.Value);
            [v1_, g1] = update_l21_spmd(v1_, x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :), weights1_, beta1.Value, Iq, ...
                dims_q, I_overlap_q, dims_overlap_q, offsetp.Value, status_q, ...
                nlevelp.Value, waveletp.Value, Ncoefs_q, temLIdxs_q, temRIdxs_q, offsetLq, offsetRq, dims_overlap_ref_q);
            g = zeros(size(x_overlap));
%             g(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :) = sigma00.Value*g0;
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
            [v2_, g2, proj_, norm_residual_check_i, norm_epsilon_check_i] = update_data_fidelity(v2_, yp, xhat_i, proj_, Ap, Atp, Gp, Wp, pUp, epsilonp, ...
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
    if ~mod(t,10)
        
        %% compute value of the priors in parallel
        spmd
            if labindex <= Qp.Value
                % compute values for the prior terms
                x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
                x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);               
                [l21_norm, nuclear_norm] = prior_overlap_spmd_cst3_weighted(x_overlap, Iq, ...
                    offset, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
                    offsetLq, offsetRq, crop_l21, crop_nuclear, w);
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
        
        % SNR
        % get xsol back from the workers
        for q = 1:Q
            xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
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
    if (param.step_flag && t>500) %rel_fval(t) < param.reweight_rel_obj)
        reweight_steps = (t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size));
        param.step_flag = 0;
    end
    
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
            (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
            t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);

        % SNR
        % get xsol back from the workers
        for q = 1:Q
            xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
        end
        sol = reshape(xsol(:),numel(xsol(:))/c,c);
        SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
        psnrh = zeros(c,1);
        for i = 1:c
            psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
        end
        SNR_average = mean(psnrh);
        
        % update weights
        spmd
            if labindex <= Qp.Value
                %x_overlap = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
                %x_overlap = zeros([dims_oq, size(xsol_q, 3)]);
                x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
                x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);

                % TO BE MODIFIED (see which element needs to be cropped, and how)
%                 [weights1_, weights0_] = update_weights_overlap(x_overlap, size(v1_), ...
%                     Iq, dims_q, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, ...
%                     Ncoefs_q, dims_overlap_ref_q, offsetLq, offsetRq, reweight_alphap);
                [weights1_, weights0_] = update_weights_overlap2_weighted(x_overlap, size(v1_), ...
                    Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, ...
                    Ncoefs_q, dims_overlap_ref_q, offsetLq, offsetRq, ...
                    reweight_alphap, crop_l21, crop_nuclear, w);    
                reweight_alphap = reweight_alpha_ffp.Value * reweight_alphap;
            else
                % compute residual image on the data nodes
                res_ = compute_residual_images(xsol(:,:,c_chunks{labindex-Qp.Value}), yp, Gp, Ap, Atp, Wp);
            end
        end
        reweight_alpha = param.reweight_alpha_ff .* reweight_alpha; % on the master node
        
        % Compute residual images (master node)
        res = zeros(size(xsol));
        for k = 1 : K
            res(:,:,c_chunks{k}) = res_{Q+k};
        end
        
        if (reweight_step_count == 0) || (reweight_step_count == 1) || (~mod(reweight_step_count,5))
            mkdir('./results/')
            save(['./results/result_HyperSARA_spmd4_new'  num2str(param.gamma) '_' num2str(reweight_step_count) '.mat'],'-v7.3','xsol','res','SNR','SNR_average');
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

% Collect distributed values (reweight_alpha, weights0_, weights1_, v0_, v1_)
v0 = cell(Q, 1);
v1 = cell(Q, 1);
weights0 = cell(Q, 1);
weights1 = cell(Q, 1);
for q = 1:Q
    v0{q} = v0_{q};
    v1{q} = v1_{q};
    weights0{q} = weights0_{q};
    weights1{q} = weights1_{q};
    xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :) = xsol_q{q};
end

% Calculate residual images
res = zeros(size(xsol));
spmd
    if labindex > Qp.Value
        res_ = compute_residual_images(xsol(:,:,c_chunks{labindex-Qp.Value}), yp, Gp, Ap, Atp, Wp);
    end
end
for k = 1 : K
    res(:,:,c_chunks{k}) = res_{Q+k};
end


norm_res = norm(res(:));
v2 = cell(K, 1);
proj = cell(K, 1);
epsilon = cell(K, 1);
for i = 1:K
    v2{i} = v2_{Q+i};
    proj{i} = proj_{Q+i};
    epsilon{i} = epsilonp{Q+i};
end

%Final log (merge this step with the computation of the residual image for
% each frequency of interest)
spmd
    if labindex <= Qp.Value
        % TO BE MODIFIED
        %x_overlap = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
        %x_overlap = zeros([dims_oq, size(xsol_q, 3)]);
        
        x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
        x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
        
        %         [l21_norm, nuclear_norm] = prior_overlap_spmd(x_overlap, Iq, ...
        %             dims_q, offset, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
        %             offsetLq, offsetRq); % d < d_sdwt2
        %         [l21_norm, nuclear_norm] = prior_overlap_spmd_cst2(x_overlap, Iq, ...
        %             dims_q, offset, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
        %             offsetLq, offsetRq, crop); % d > d_sdwt2
        [l21_norm, nuclear_norm] = prior_overlap_spmd_cst3_weighted(x_overlap, Iq, ...
            offset, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
            offsetLq, offsetRq, crop_l21, crop_nuclear, w);
    end
end

l21 = 0;
nuclear = 0;
for q = 1:Q
    l21 = l21 + l21_norm{q};
    nuclear = nuclear + nuclear_norm{q};
end

% SNR (only on the master node this time)
sol = reshape(xsol(:),numel(xsol(:))/c,c);
SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
psnrh = zeros(c,1);
for i = 1:c
    psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
end
SNR_average = mean(psnrh);

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
