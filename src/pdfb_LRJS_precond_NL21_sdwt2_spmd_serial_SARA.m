function [xsol,v0,v1,v2,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = pdfb_LRJS_precond_NL21_sdwt2_spmd_serial_SARA(y, epsilon, A, At, pU, G, W, param, X0, K, wavelet, nlevel, c_chunks, c)

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

% total number of workers (Q: facets workers, K: data workers)
numworkers = 2 + K;
cirrus_cluster = parcluster('local');
cirrus_cluster.NumWorkers = numworkers;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
if cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads > ncores
    exit(1);
end

parpool(cirrus_cluster, numworkers);

% instantiate Psi, Psit
[Psi, Psit] = op_sp_wlt_basis(wavelet, nlevel, M, N);
spmd
    if labindex == 2
        [Psi_, Psit_] = op_sp_wlt_basis(wavelet, nlevel, M, N);
    end
end

% Initializations.
if isfield(param,'init_xsol')
    xsol = param.init_xsol;
    fprintf('xsol uploaded \n\n')
else
    xsol = zeros(M,N,c);
    fprintf('xsol NOT uploaded \n\n')
end
g = zeros(size(xsol));

% Primal / prior nodes (l21/nuclear norm dual variables)
v0_ = Composite();
weights0_ = Composite();
v1_ = Composite();
weights1_ = Composite();
if isfield(param,'init_v0') || isfield(param,'init_v1')
    v0_{1} = param.init_v0{1};
    v1_{2} = param.init_v1{1};
    weights0_{1} = param.init_weights0{1};
    weights1_{2} = param.init_weights1{1};
    s = size(param.init_v1{1}, 1);
else
    spmd
        if labindex == 1
            v0_ = zeros(M*N, c);
            weights0_ = ones(min(M*N, c), 1);
        elseif labindex == 2
            [v1_, weights1_, s] = initialize_l21_serial(xsol, Psit_, 'zpd', nlevel);
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

% to be cleansed later on (change the format of the input data?)
% see warm-restart in detail for this step...
Ap = Composite();
Atp = Composite();
Gp = Composite();
yp = Composite();
pUp = Composite();
Wp = Composite();
epsilonp = Composite();
norm_res = Composite();
xhat_i = Composite();
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
    yp{2+k} = y{k};
    y{k} = [];
    Ap{2+k} = A;
    Atp{2+k} = At;
    Gp{2+k} = G{k};
    G{k} = [];
    Wp{2+k} = W{k};
    W{k} = [];
    pUp{2+k} = pU{k};
    pU{k} = [];
    epsilonp{2+k} = epsilon{k};
    norm_res{2+k} = norm_res_tmp;
end

clear norm_res_tmp epsilon pU W G y

v2_ = Composite();
t_block = Composite();
proj_ = Composite();
if isfield(param,'init_v2') % assume all the other related elements are also available in this case
    for k = 1:K
        v2_{2+k} = param.init_v2{k};
        proj_{2+k} = param.init_proj{k};
        t_block{2+k} = param.init_t_block{k};
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
        v2_{2+k} = v2_tmp;
        proj_{2+k} = proj_tmp;
        t_block{2+k} = t_block_;
    end
end

clear proj_tmp v2_tmp norm_res_tmp t_block_ G y

reweight_last_step_iter = 0;
reweight_step_count = 0;
rw_counts = 1;

%% Reweighting parameters

reweight_alpha = param.reweight_alpha;
reweight_steps = param.reweight_steps;

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
    
    %% Update primal variable
    prev_xsol = xsol;
    xsol = max(real(xsol - g), 0);
    xhat = 2*xsol - prev_xsol;

    for i = 1:K
       xhat_i{2+i} = xhat(:, :, c_chunks{i}); 
    end
    
    %% Update dual variables
    spmd
        if labindex == 1 % nuclear prior node
            tic;
            [v0_, g0_] = run_par_nuclear(v0_, xhat, weights0_, beta0.Value, sigma00.Value);
            tw = toc;  
        elseif labindex == 2 % l21 prior node (full SARA prior)
            tic
            [v1_, g1_] = run_par_l21(v1_, Psit_, Psi_, xhat, weights1_, beta1.Value, sigma11.Value);  
            tw = toc;
        else % data nodes, 3:K+2 (labindex > 2)
            tic
            [v2_, g2_, proj_, norm_residual_check_i, norm_epsilon_check_i] = update_data_fidelity(v2_, yp, xhat_i, proj_, Ap, Atp, Gp, Wp, pUp, epsilonp, ...
                elipse_proj_max_iter.Value, elipse_proj_min_iter.Value, elipse_proj_eps.Value, sigma22.Value);
            tw = toc;
        end
    end
    
    g = g0_{1} + g1_{2};
    for i = 1:K
        g(:,:,c_chunks{i}) = g(:,:,c_chunks{i}) + g2_{2+i};
    end
    
    %% Relative change of objective function
    rel_x = norm(prev_xsol(:) - xsol(:));
    norm_x = norm(xsol(:));
    rel_fval(t) = rel_x/norm_x;
    end_iter(t) = toc(start_iter);
    fprintf('Iter = %i, Time = %e\n',t,end_iter(t));
    
%     t_nuclear = tw{1};
%     t_l21 = tw{2};
%     t_data = 0;
%     for k = 1:K
%        t_data = max(t_data, tw{2+k}); 
%     end
    
    %% Display
    if ~mod(t,10)
        
        % [P.-A.]
        %% compute value of the priors
        % nuclear norm
        nuclear = nuclear_norm(xsol);
        
        % l21 norm
        l21 = l21_norm_sara(xsol, Psit, s{2});
        
        
        % retrieve value of the monitoring variables (residual norms)
        norm_epsilon_check = 0;
        norm_residual_check = 0;
        for i = 3:K+2
            norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{i};
            norm_residual_check = norm_residual_check + norm_residual_check_i{i};
        end
        norm_epsilon_check = sqrt(norm_epsilon_check);
        norm_residual_check = sqrt(norm_residual_check);
        
        % SNR
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
            if labindex > 2
                [epsilonp, t_block] = update_epsilon(epsilonp, t, t_block, rel_fval(t), norm_res, ...
                    adapt_eps_tol_in.Value, adapt_eps_tol_out.Value, adapt_eps_steps.Value, adapt_eps_rel_obj.Value, ...
                    adapt_eps_change_percentage.Value);
            end
        end
    end
    
    %% Reweighting (in parallel)
    if (param.step_flag && rel_fval(t) < param.reweight_rel_obj)
        reweight_steps = (t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size));
        param.step_flag = 0;
    end
    
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
            (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
            t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);
        
        spmd
            if labindex == 1
                weights0_ = update_weights_nuclear_serial(xsol, reweight_alpha);
            elseif labindex == 2
                weights1_ = update_weights_l21_serial(xsol, Psit_, weights1_, reweight_alpha);
            end
        end
        reweight_alpha = param.reweight_alpha_ff .* reweight_alpha;
        
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
v0 = v0_{1};
v0_{1} = [];
v1 = v1_{2};
v1_{2} = [];
weights0 = weights0_{1};
weights1 = weights1_{2};

% to be completely modified (within spmd function?)
% Calculate residual images
res = zeros(size(xsol));
spmd
    if labindex > 2
        res_ = compute_residual_images(xsol(:,:,c_chunks{labindex-2}), yp, Gp, Ap, Atp, Wp);
    end
end
for k = 1 : K
    res(:,:,c_chunks{k}) = res_{2+k};
    res_{2+k} = [];
end

norm_res = norm(res(:));
v2 = cell(K, 1);
proj = cell(K, 1);
epsilon = cell(K, 1);
for i = 1:K
    v2{i} = v2_{2+i};
    v2_{2+i} = [];
    proj{i} = proj_{2+i};
    proj_{2+i} = [];
    epsilon{i} = epsilonp{2+i};
    epsilonp{2+i} = [];
end

%Final log (merge this step with the computation of the residual image for
% each frequency of interest)
% nuclear norm
nuclear = nuclear_norm(xsol);
% l21 norm
l21 = l21_norm_sara(xsol, Psit, s{2});

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
