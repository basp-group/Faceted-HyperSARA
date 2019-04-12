function [xsol,v0,v1,v2,weights0,weights1,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = pdfb_LRJS_precond_NL21_sdwt2_spmd4(y, epsilon, A, At, pU, G, W, param, X0, Qx, Qy, K, wavelet, L, nlevel, c_chunks, c)

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

% oversampling vectorized data length
No = size(W{1}{1}{1}, 1);

% number of pixels
[M, N] = size(At(zeros(No, 1)));
x0 = reshape(X0, [M, N, c]);

% [P.-A.]
% define spatial facets (no overlap)
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
% [P.-A.] instantiate auxiliary variables for sdwt2
[~, ~, ~, dims_overlap_ref, I_overlap, dims_overlap, ...
    ~, ~, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices([M, N], I, dims, nlevel, wavelet, L);
id_dirac = find(ismember(wavelet, 'self'), 1);
dirac_present = ~isempty(id_dirac);

% [P.-A.]
% total number of workers (Q: facets workers, K: data workers)
numworkers = Q + K;
cirrus_cluster = parcluster('local');
cirrus_cluster.NumWorkers = numworkers;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
if cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads > ncores
    exit(1);
end

% [P.-A.] /!\ only simple indexing allowed into Composite objects from the master
parpool(cirrus_cluster, numworkers); % override default preference

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
if isfield(param,'init_xsol')
    xsol = param.init_xsol;
    fprintf('xsol uploaded \n\n')
else
    xsol = zeros(M,N,c);
    fprintf('xsol NOT uploaded \n\n')
end

% Prior/primal nodes
% l21 / nuclear norm dual variables
v0_ = Composite();
weights0_ = Composite();
v1_ = Composite();
weights1_ = Composite(); % idem, can be created in parallel with spmd
for q = 1:Q
    p = prod(Ncoefs{q}, 2);
    if dirac_present
        sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) + prod(dims(q,:));
    else
        sz = 3*sum(p) - 2*sum(p(nlevel+1:nlevel+1:end));
    end
    
    v0_{q} = zeros(prod(dims_overlap_ref(q,:)), c);
    weights0_{q} = zeros(min(prod(dims_overlap_ref(q, :)), c), 1);
    v1_{q} = zeros(sz, c);
    weights1_{q} = zeros(sz, c);
end
clear sz

%% Data node parameters: to be completely modified
% data nodes
% to be done properly with an spmd

A = afclean(A);
At = afclean(At);

% see if these constants are properly updated in practice
elipse_proj_max_iter = parallel.pool.Constant(param.elipse_proj_max_iter);
elipse_proj_min_iter = parallel.pool.Constant(param.elipse_proj_min_iter);
elipse_proj_eps = parallel.pool.Constant(param.elipse_proj_eps);
adapt_eps_tol_in = parallel.pool.Constant(param.adapt_eps_tol_in);
adapt_eps_steps = parallel.pool.Constant(param.adapt_eps_steps);
adapt_eps_rel_obj = parallel.pool.Constant(param.adapt_eps_rel_obj);
adapt_eps_change_percentage = parallel.pool.Constant(param.adapt_eps_change_percentage);

% count_eps_update_down, count_eps_update_up

Ap = Composite();
Atp = Composite();
x_hat_i = Composite();
v2_ = Composite();
t_block = Composite();
Gp = Composite();
yp = Composite();
pUp = Composite();
Wp = Composite();
epsilonp = Composite();
norm_res = Composite();
proj = Composite();
for k = 1:K
    v2_tmp = cell(length(c_chunks{k}), 1);
    t_block_ = cell(length(c_chunks{k}), 1);
    norm_res_tmp = cell(length(c_chunks{k}), 1);
    proj_tmp = cell(length(c_chunks{k}), 1);
    for i = 1:length(c_chunks{k})
        v2_tmp{i} = cell(length(G{k}{i}),1);
        t_block_{i} = cell(length(G{k}{i}),1);
        norm_res_tmp{i} = cell(length(G{k}{i}),1);
        proj_tmp{i} = cell(length(G{k}{i}),1);
        for j = 1 : length(G{k}{i})
            v2_tmp{i}{j} = zeros(length(y{k}{i}{j}) ,1);
            t_block_{i}{j} = 0;
            norm_res_tmp{i}{j} = norm(y{k}{i}{j});
            proj_tmp{i}{j} = zeros(size(y{k}{i}{j}));
        end
    end
    v2_{Q+k} = v2_tmp;
    yp{Q+k} = y{k};
    x_hat_i{Q+k} = zeros(M, N, length(c_chunks{k}));
    Ap{Q+k} = A;
    Atp{Q+k} = At;
    Gp{Q+k} = G{k};
    Wp{Q+k} = W{k};
    pUp{Q+k} = pU{k};
    epsilonp{Q+k} = epsilon{k};
    norm_res{Q+k} = norm_res_tmp;
    proj{Q+k} = proj_tmp;
end

clear proj_tmp v2_tmp norm_res_tmp t_block_

t_start = 1;
reweight_last_step_iter = 0;
reweight_step_count = 0;
rw_counts = 1;

% To be defined
% t, t_block, rel_fval, norm_res

%% Reweighting parameters

% [P.-A.]
reweight_alpha = param.reweight_alpha;
reweight_alphap = Composite();
for q = 1:Q
    reweight_alphap{q} = reweight_alpha;
end
reweight_alpha_ff = parallel.pool.Constant(param.reweight_alpha_ff);
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

%Step size primal
tau = 0.99/(sigma0*param.nu0 + sigma1*param.nu1 + sigma2*param.nu2);

sigma00 = parallel.pool.Constant(tau*sigma0);
sigma11 = parallel.pool.Constant(tau*sigma1);
sigma22 = parallel.pool.Constant(tau*sigma2);

flag = 0;

% [P.-A.]
beta0 = parallel.pool.Constant(param.gamma0/sigma0); % needed only on the "prior" workers
beta1 = parallel.pool.Constant(param.gamma/sigma1);

% Main loop. Sequential.
% util_create_pool(param.num_workers);
% maxNumCompThreads(param.num_workers);

start_loop = tic;

for t = t_start : param.max_iter
    
    %fprintf('Iter %i\n',t);
    start_iter = tic;
    
    spmd
        if labindex <= Q
            % primal/prior nodes (1:Q)
            % update primal variable
            [xsol_q, xhat_q, rel_x_q, norm_x_q] = run_par_primal(xsol_q, g_q);
            
            % send xhat_q (communication towards the data nodes)
            for i = 1:K
                labSend(xhat_q(:,:,c_chunksp.Value{i}), Q+i); % problem indexing here!
            end
            
            % update ghost cells (versions of xhat with overlap)
            % overlap_q = dims_overlap_ref_q - dims_q;
            x_overlap = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
            x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xhat_q;
            x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp, Qxp); % problem index in l. 80 (position 2...) on worker 2 (to be investigated further)
            
            % update dual variables (nuclear, l21)
            [v0_, g0] = run_par_nuclear_spmd(v0_, x_overlap, weights0_, beta0.Value);
            [v1_, g1] = run_par_l21_spmd(v1_, x_overlap, weights1_, beta1.Value, Iq, ...
                dims_q, I_overlap_q, dims_overlap_q, offsetp.Value, status_q, ...
                nlevelp.Value, waveletp.Value, Ncoefs_q, temLIdxs_q, temRIdxs_q, offsetLq, offsetRq, dims_overlap_ref_q);
            g1 = comm2d_reduce_wideband(g1, overlap, Qyp, Qxp, Kp);
            
            % compute g_ for the final update term
            g_ = sigma00.Value*g0(overlap(1)+1:end, overlap(2)+1:end, :) + ...
                sigma11.Value*g1(overlap(1)+1:end, overlap(2)+1:end, :);
            
            % retrieve portions of g2 from the data nodes
            for i = 1:Kp.Value
                g_(:,:,c_chunksp.Value{i}) = g_(:,:,c_chunksp.Value{i}) + labReceive(Qp.Value+i);
            end
        else
            % data nodes (Q+1:Q+K) (no data blocking, just frequency for the moment
            % retrieve xhat_i from the prior/primal nodes
            for q = 1:Qp.Value
                xhat_i(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2), :) = ...
                    labReceive(q);
            end
            [v2_, g2, norm_residual_check_i, norm_epsilon_check_i] = run_par_data_fidelity(v2_, yp, xhat_i, proj, Ap, Atp, Gp, Wp, pUp, epsilonp, ...
                elipse_proj_max_iter.Value, elipse_proj_min_iter.Value, elipse_proj_eps.Value, sigma22.Value); % error somewhere inside the code
            
            % send portions of g2 to the prior/primal nodes
            for q = 1:Qp.Value
                labSend(g2(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2), :), q); % to be done inside the function?
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
        
        % [P.-A.]
        %% compute value of the priors in parallel (include weights*?)
        spmd
            if labindex <= Q
                % compute values for the prior terms
                %x_overlap = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
                %see if this is necessary or not (to be further investigated)
                x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
                x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp, Qxp);

                [l21_norm, nuclear_norm] = prior_value_spmd(x_overlap, overlap, Iq, ...
                    dims_q, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
                    offsetLq, offsetRq);
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
        
%         figure(1),
%         subplot(2,2,1);
%         imagesc(log10(max(flip(x0(:,:,1)),0))); hold on; colorbar; axis image; axis off; colormap(cubehelix); caxis([-3.5, 0]);
%         subplot(2,2,2);
%         imagesc(log10(max(flip(xsol(:,:,1)),0))); hold on; colorbar; axis image; axis off; colormap(cubehelix); caxis([-3.5, 0]);
%         subplot(2,2,3);
%         imagesc(log10(max(flip(x0(:,:,end)),0))); hold on; colorbar; axis image; axis off; colormap(cubehelix); caxis([-3.5, 0]);
%         subplot(2,2,4);
%         imagesc(log10(max(flip(xsol(:,:,end)),0))); hold on; colorbar; axis image; axis off; colormap(cubehelix); caxis([-3.5, 0]);
%         pause(0.1)
        
        fitswrite(x0, './x0.fits');
        fitswrite(xsol, './xsol.fits');
    end
    
    %% Global stopping criteria
    %     prod(prod(residual_check < param.adapt_eps_tol_out*epsilon_check)) && prod(prod(residual_check > param.adapt_eps_tol_in*epsilon_check))
    if t>1 && rel_fval(t) < param.rel_obj && reweight_step_count > param.total_reweights && ...
            (norm_residual_check <= param.adapt_eps_tol_out*norm_epsilon_check)
        flag = 1;
        break;
    end
    
    %% Update epsilons (to be completely modified) (to be done inside an spmd, in parallel of the reweighting procedure)
    % write function update_epsilon
    if param.use_adapt_eps && t > param.adapt_eps_start
        spmd
            if labindex > Q
                epsilonp = update_epsilon(epsilonp, t, t_block, rel_fval(t), norm_res, ...
                    adapt_eps_tol_in.Value, adapt_eps_steps.Value, adapt_eps_rel_obj.Value, ...
                    adapt_eps_change_percentage.Value); % use some variable sonly defined on the master node... see if the communication works as intendd in this case...
            end
        end
    end
    
    %% Reweighting (to be partly modified)
    if (param.step_flag && rel_fval(t) < param.reweight_rel_obj) % && (norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check))
        reweight_steps = [t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size)];
        param.step_flag = 0;
    end
    
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
            (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
            t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);
        
        % [P.-A.]
        spmd
            if labindex <= Qp.Value
                x_overlap = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
                x_overlap(overlap(1)+1:end, overlap(2)+1:end) = xsol_q;
                x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp, Qxp);

                % nuclear norm
                sol = reshape(xsol_q, [size(xsol_q, 1)*size(xsol_q, 2), size(xsol_q, 3)]);
                [~,S00,~] = svd(sol,'econ');
                d_val0 = abs(diag(S00));
                weights0_ = reweight_alpha ./ (reweight_alpha + d_val0);

                % l21 norm
                zerosNum = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding
                x_ = zeros([zerosNum, size(x_overlap, 3)]);
                x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap;
                w = zeros(size(v1_));
                for l = 1 : size(x_, 3)
                    w(:, l) = sdwt2_sara(x_(:, :, l), Iq, dims_q, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q);
                end
                d_val1 = sqrt(sum(abs((w)).^2,2));
                weights1_ = reweight_alphap ./ (reweight_alphap + d_val1);
                reweight_alphap = reweight_alpha_ff .* reweight_alphap;
%             else
%                 % compute residual image on the data nodes
%                 res_ = compute_residual_images(x, y, G, A, At); % to be possibly removed (not) used here...
            end
        end
        reweight_alpha = reweight_alpha_ff .* reweight_alpha; % on the master node
        %-- end modifications
        
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
end_loop = toc(start_loop)

% [P.-A.] collect distributed values (reweight_alpha,weights0, weights1)
% do sth for the weights... (to be done later on)
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

% to be completely modified (within spmd function?)
% Calculate residual images
res = zeros(size(xsol));
for k = 1 : K
    for i = 1 : length(c_chunks{k})
        Fx = A(xsol(:,:,c_chunks{k}(i)));
        g2 = zeros(No,1);
        for j = 1 : length(G{k}{i})
            res_f = y{k}{i}{j} - G{k}{i}{j} * Fx(W{k}{i}{j});
            u2 = G{k}{i}{j}' * res_f;
            g2(W{k}{i}{j}) = g2(W{k}{i}{j}) + u2; % no overlap between different groups? (content of W...)
        end
        res(:,:,c_chunks{k}(i)) = real(At(g2)); % residual images
    end
end

norm_res = norm(res(:));
v2 = cell(K, 1);
for i = 1:K
    v2{i} = v2_{Q+i};
end

%Final log (merge this step with the computation of the residual image for
% each frequency of interest)
spmd
    x_overlap = zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
    x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xsol_q;
    x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp, Qxp);
    
    [l21_norm, nuclear_norm] = prior_value_spmd(x_overlap, overlap, Iq, ...
        dims_q, offset, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, dims_overlap_ref_q, ...
        offsetLq, offsetRq);
end

l21 = 0;
nuclear = 0;
for q = 1:Q
    l21 = l21 + l21_norm{q};
    nuclear = nuclear + nuclear_norm{q};
end

% SNR (only on the master node this time? to be seen)
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
        fprintf('epsilon = %e, residual = %e\n', norm(epsilon_check),norm(residual_check));
        fprintf('SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf('Iter %i\n',t);
        fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
        fprintf('epsilon = %e, residual = %e\n', norm(epsilon_check),norm(residual_check));
        fprintf('SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
    end
end

end
