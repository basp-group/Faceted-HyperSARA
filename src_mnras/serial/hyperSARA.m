function [xsol,param,epsilon,t,rel_val,nuclear,l21,norm_res_out,res,end_iter,SNR,SNR_average] = ...
    hyperSARA(y, epsilon, A, At, pU, G, W, param, X0, K, wavelet, nlevel, c_chunks, c, init_file_name, name)

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
norm_epsilon_check = Inf;
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

%Initializations.
init_flag = isfile(init_file_name)
if init_flag
    init_m = matfile(init_file_name)
end

if init_flag
    xsol = init_m.xsol;
    % xhat = init_m.xhat;
    param = init_m.param;
    epsilon = init_m.epsilon;
    fprintf('xsol, param and epsilon uploaded \n\n')
    %! display norm of interest
%     norm_eps = 0;
%     for k = 1:K
%         for l = 1:numel(epsilon{k})
%             for b = 1:numel(epsilon{k}{l})
%                 norm_eps = norm_eps + norm(epsilon{k}{l}{b}, 'fro');
%             end
%         end
%     end
%     fprintf('|xsol| = %e, |epsilon| = %e \n\n', norm(xsol(:)), norm_eps)
%     param
else
    xsol = zeros(M,N,c);
    fprintf('xsol initialized \n\n')
end

% Primal / prior nodes (l21/nuclear norm dual variables)
v0_ = Composite();
weights0_ = Composite();
v1_ = Composite();
s_ = Composite();
weights1_ = Composite();
if init_flag
    v0_{1} = init_m.v0;
    v1_{2} = init_m.v1;
    s_{2} = size(v1_{2});
    weights0_{1} = init_m.weights0;
    weights1_{2} = init_m.weights1;
    s = size(init_m.v1, 1);
    fprintf('v0, v1, weigths0, weights1 uploaded \n\n')
    %!
%     fprintf('|v0| = %e, |v1| = %e, |weights0| = %e, |weights1| = %e \n\n', norm(init_m.v0,'fro'), norm(init_m.v1,'fro'), ...
%         norm(init_m.weights0,'fro'), norm(init_m.weights1,'fro'))
else
    spmd
        if labindex == 1
            v0_ = zeros(M*N, c);
            weights0_ = ones(min(M*N, c), 1);
        elseif labindex == 2
            [v1_, weights1_, s_] = initialize_l21_serial(xsol, Psit_, 'zpd', nlevel);
        end
    end
    s = s_{2};
    fprintf('v0, v1, weigths0, weights1 initialized \n\n')
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
adapt_eps_rel_var = parallel.pool.Constant(param.adapt_eps_rel_var);
adapt_eps_change_percentage = parallel.pool.Constant(param.adapt_eps_change_percentage);

Ap = Composite();
Atp = Composite();
xhat_i = Composite();
Gp = Composite();
yp = Composite();
pUp = Composite();
Wp = Composite();
epsilonp = Composite();
% [09/10/2019] fixing bug in initialization of norm_res (warm-restart)
norm_res = Composite();
if init_flag
    for k = 1:K
        norm_res(2+k) = init_m.norm_res(k,1);
    end
    fprintf('norm_res uploaded \n\n')
    
    %!
%     nres = 0;
%     for k = 1:K
%         e = norm_res{2+k};
%         for l = 1:numel(e)
%             for b = 1:numel(e{l})
%                 nres = nres + norm(e{l}{b})^2;
%             end
%         end
%     end
%     nres = sqrt(nres);
%     fprintf('|norm_res| = %e \n\n', nres);
else
    for k = 1:K
        norm_res_tmp = cell(length(c_chunks{k}), 1);
        for i = 1:length(c_chunks{k})
            norm_res_tmp{i} = cell(length(y{k}{i}),1);
            for j = 1 : length(y{k}{i})
                norm_res_tmp{i}{j} = norm(y{k}{i}{j});
            end
        end
        norm_res{2+k} = norm_res_tmp;
    end
    clear norm_res_tmp
    fprintf('norm_res initialized \n\n')
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
end
clear epsilon pU W G y

v2_ = Composite();
t_block = Composite();
proj_ = Composite();
if init_flag
    for k = 1:K
        v2_(2+k) = init_m.v2(k,1);
        proj_(2+k) = init_m.proj(k,1);
        t_block(2+k) = init_m.t_block(k,1);
    end
    fprintf('v2, proj, t_block uploaded \n\n')
    %!
%     nv2 = 0;
%     nproj = 0;
%     ntblock = 0;
%     for k = 1:K
%         e = init_m.v2(k,1);
%         e1 = init_m.proj(k,1);
%         e2 = init_m.t_block(k,1);
%         for l = 1:numel(e{1})
%             for b = 1:numel(e{1}{l})
%                 nv2 = nv2 + norm(e{1}{l}{b}, 'fro');
%                 nproj = nproj + norm(e1{1}{l}{b}, 'fro');
%                 ntblock = ntblock + norm(e2{1}{l}{b}, 'fro');
%             end
%         end
%     end
%     fprintf('|v2| = %e, |proj| = %e, |t_block| = %e \n\n', nv2, nproj, ntblock);
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
    clear proj_tmp v2_tmp t_block_
    fprintf('v2, proj, t_block initialized \n\n')
end

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
    reweight_last_step_iter = 0;
    fprintf('reweight_last_iter_step initialized \n\n')
end
rw_counts = 1;

%% Reweighting parameters
% sig_bar = param.sig_bar;
% sig = param.sig;
reweight_alpha = param.reweight_alpha;
reweight_alphap = Composite();
for q = 1:2
    reweight_alphap{q} = reweight_alpha;
end
reweight_alpha_ffp = parallel.pool.Constant(param.reweight_alpha_ff);
reweight_steps = param.reweight_steps;

if init_flag
    g = init_m.g;
    fprintf('g uploaded \n\n')
    %!
    % fprintf('|g| = %e \n\n', norm(g(:)))
else
    g = zeros(size(xsol));
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

if isfield(param, 'init_t_start')
    t_start = param.init_t_start;
    fprintf('t_start uploaded \n\n')
else
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
    rel_val = zeros(param.max_iter, 1);
    end_iter = zeros(param.max_iter, 1);
    t_master = zeros(param.max_iter, 1);
    t_l21 = zeros(param.max_iter, 1);
    t_nuclear = zeros(param.max_iter, 1);
    t_data = zeros(param.max_iter, 1);
    fprintf('rel_val, end_iter, t_master, t_l21, t_nuclear and t_data initialized \n\n')
end

%%
%! check warm-start worked as expected
if init_flag
    spmd
        if labindex > 2
            [norm_residual_check_i, norm_epsilon_check_i] = sanity_check(epsilonp, norm_res);
        end
    end  
    norm_epsilon_check = 0;
    norm_residual_check = 0;
    for i = 3:K+2
        norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{i};
        norm_residual_check = norm_residual_check + norm_residual_check_i{i};
    end
    norm_epsilon_check = sqrt(norm_epsilon_check);
    norm_residual_check = sqrt(norm_residual_check);
    
    % nuclear norm
    nuclear = nuclear_norm(xsol);
    
    % l21 norm
    l21 = compute_sara_prior(xsol, Psit, s);
    
    % SNR
    sol = reshape(xsol(:),numel(xsol(:))/c,c);
    SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
    psnrh = zeros(c,1);
    for i = 1:c
        psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
    end
    SNR_average = mean(psnrh);
    
    % Log
    if (param.verbose >= 1)
        fprintf('Iter %i\n',t_start-1);
        fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t_start-1));
        fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
        fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
    end
end
%%

start_loop = tic;
fprintf('START THE LOOP MNRAS ver \n\n')
param.max_iter = 100000;

for t = t_start : param.max_iter
    
    %fprintf('Iter %i\n',t);
    start_iter = tic;
    
    %% Update primal variable
    tw = tic;
    prev_xsol = xsol;
    xsol = max(real(xsol - g), 0);
    xhat = 2*xsol - prev_xsol;
    t_master(t) = toc(tw);

    for i = 1:K
       xhat_i{2+i} = xhat(:, :, c_chunks{i}); 
    end
    
    %% Update dual variables
    spmd
        if labindex == 1 % nuclear prior node
            tw = tic;
            [v0_, g0_] = update_dual_nuclear_serial(v0_, xhat, weights0_, beta0.Value, sigma00.Value);
            t_op = toc(tw);  
        elseif labindex == 2 % l21 prior node (full SARA prior)
            tw = tic;
            [v1_, g1_] = update_dual_l21_serial(v1_, Psit_, Psi_, xhat, weights1_, beta1.Value, sigma11.Value);  
            t_op = toc(tw); 
        else % data nodes, 3:K+2 (labindex > 2)
            tw = tic;
            [v2_, g2_, proj_, norm_res, norm_residual_check_i, norm_epsilon_check_i] = update_dual_fidelity(v2_, yp, xhat_i, proj_, Ap, Atp, Gp, Wp, pUp, epsilonp, ...
                elipse_proj_max_iter.Value, elipse_proj_min_iter.Value, elipse_proj_eps.Value, sigma22.Value);
            t_op = toc(tw); 
        end
    end
    
    g = g0_{1} + g1_{2};
    for i = 1:K
        g(:,:,c_chunks{i}) = g(:,:,c_chunks{i}) + g2_{2+i};
    end
    
    %% Relative change of objective function
    rel_x = norm(prev_xsol(:) - xsol(:));
    norm_x = norm(xsol(:));
    rel_val(t) = rel_x/norm_x;
    end_iter(t) = toc(start_iter);
    
    % retrieve update time (average for data processes)
    t_nuclear(t) = t_op{1};
    t_l21(t) = t_op{2};
    t_data(t) = 0; % just in case
    for i = 3:K+2
        t_data(t) = t_data(t) + t_op{i};
    end
    t_data(t) = t_data(t)/K;
    fprintf('Iter = %i, Time = %e, t_master= %e, t_l21 = %e, t_nuclear = %e, t_data = %e\n',t,end_iter(t),t_master(t),t_l21(t),t_nuclear(t),t_data(t));
    
    %% Retrieve value of the monitoring variables (residual norms)
    norm_epsilon_check = 0;
    norm_residual_check = 0;
    for i = 3:K+2
        norm_epsilon_check = norm_epsilon_check + norm_epsilon_check_i{i};
        norm_residual_check = norm_residual_check + norm_residual_check_i{i};
    end
    norm_epsilon_check = sqrt(norm_epsilon_check);
    norm_residual_check = sqrt(norm_residual_check);
    
%     t_nuclear = tw{1};
%     t_l21 = tw{2};
%     t_data = 0;
%     for k = 1:K
%        t_data = max(t_data, tw{2+k}); 
%     end
    
    %% Display
    if ~mod(t,100)
        
        % [P.-A.]
        %% compute value of the priors
        % nuclear norm
        nuclear = nuclear_norm(xsol);
        
        % l21 norm
        l21 = compute_sara_prior(xsol, Psit, s);
        
        % SNR
        sol = reshape(xsol(:),numel(xsol(:))/c,c);
        SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
        psnrh = zeros(c,1);
        for i = 1:c
            psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
        end
        SNR_average = mean(psnrh);
        
        % Log
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
            fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
            fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
        end 
    end
    
    %% Global stopping criteria
    if t>1 && rel_val(t) < param.rel_var && reweight_step_count > param.total_reweights && ...
            (norm_residual_check <= param.adapt_eps_tol_out*norm_epsilon_check)
    % if ((t>1) && (reweight_step_count >= param.total_reweights)) && ((rel_val(t) < param.rel_var && ...
    %     (norm(residual_check) < param.adapt_eps_tol_out*norm(epsilon_check))) || ...
    %     (t - reweight_last_step_iter >= param.ppd_max_iter))
        flag = 1;
        break;
    end
    
    %% Update epsilons (in parallel)
    if param.use_adapt_eps && (t > param.adapt_eps_start) && (rel_val(t) < param.adapt_eps_rel_var)
        spmd
            if labindex > 2
                [epsilonp, t_block] = update_epsilon(epsilonp, t, t_block, rel_val(t), norm_res, ...
                    adapt_eps_tol_in.Value, adapt_eps_tol_out.Value, adapt_eps_steps.Value, adapt_eps_rel_var.Value, ...
                    adapt_eps_change_percentage.Value);
            end
        end
    end
    
    %% Reweighting (in parallel)
    if (param.step_flag && t > 500) %rel_val(t) < param.reweight_rel_var)
        reweight_steps = (t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size));
        param.step_flag = 0;
    end
        
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
            (param.use_reweight_eps && rel_val(t) < param.reweight_rel_var && ...
            norm_residual_check <= param.adapt_eps_tol_out*norm_epsilon_check && ...
            t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr) || ...
            (t - reweight_last_step_iter > 3000)
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);

        % SNR (only on the master node)
        sol = reshape(xsol(:),numel(xsol(:))/c,c);
        SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
        psnrh = zeros(c,1);
        for i = 1:c
            psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
        end
        SNR_average = mean(psnrh);
        
        spmd
            if labindex == 1
                weights0_ = update_weights_nuclear_serial(xsol, reweight_alpha);
            elseif labindex == 2
                weights1_ = update_weights_l21_serial(xsol, Psit_, weights1_, reweight_alpha);
            else % > 2
                % compute residual image
                res_ = compute_residual_images(xsol(:,:,c_chunks{labindex-2}), yp, Gp, Ap, Atp, Wp);
            end
        end
        reweight_alpha = param.reweight_alpha_ff .* reweight_alpha;
        param.reweight_alpha = reweight_alpha;
        param.init_reweight_step_count = reweight_step_count+1;
        param.init_reweight_last_iter_step = t;
        param.init_t_start = t+1; %! should  be t+1 ?

        if (reweight_step_count >= param.total_reweights)
            % param.reweight_max_reweight_itr = t+1;
            fprintf('\n\n No more reweights \n\n');
            break;
        end  
        
        if (reweight_step_count == 0) || (reweight_step_count == 1) || (~mod(reweight_step_count,5))
            m = matfile([name, '_', ...
              num2str(param.ind) '_' num2str(param.gamma) '_' num2str(reweight_step_count) '.mat'], ...
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
            m.v0 = v0_{1};
            m.v1 = v1_{2};
            m.weights0 = weights0_{1};
            m.weights1 = weights1_{2};
            % Retrieve variables from workers
            % data nodes
            for k = 1:K
                m.res(:,:,c_chunks{k}) = res_{2+k};
                res_{2+k} = [];
                m.v2(k,1) = v2_(2+k);
                m.proj(k,1) = proj_(2+k);
                m.t_block(k,1) = t_block(2+k);
                m.epsilon(k,1) = epsilonp(2+k);
                m.norm_res(k,1) = norm_res(2+k);
            end
            m.SNR = SNR;
            m.SNR_average = SNR_average;
            m.end_iter = end_iter;
            m.t_master = t_master;
            m.t_l21 = t_l21;
            m.t_nuclear = t_nuclear;
            m.t_data = t_data;
            m.rel_val = rel_val;
            clear m
            
            % nuclear norm
            nuclear = nuclear_norm(xsol);
            
            % l21 norm
            l21 = compute_sara_prior(xsol, Psit, s);
            
            % SNR
            sol = reshape(xsol(:),numel(xsol(:))/c,c);
            SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
            psnrh = zeros(c,1);
            for i = 1:c
                psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
            end
            SNR_average = mean(psnrh);
            
            % Log
            if (param.verbose >= 1)
                fprintf('Backup iter: %i\n',t);
                fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
                fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
                fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
            end
        end
        
        %%
        %! debugging
%         norm_eps = 0;
%         for k = 1:K
%             e = epsilonp{2+k};
%             for l = 1:numel(e)
%                 for b = 1:numel(e{l})
%                     norm_eps = norm_eps + norm(e{l}{b}, 'fro');
%                 end
%             end
%         end
%         fprintf('|xsol| = %e, |epsilon| = %e \n\n', norm(xsol(:)), norm_eps)
%         fprintf('|v0| = %e, |v1| = %e, |weights0| = %e, |weights1| = %e \n\n', norm(m.v0, 'fro'), norm(m.v1, 'fro'), ...
%         norm(m.weights0,'fro'), norm(m.weights1,'fro'))
%         nres = 0;
%         for k = 1:K
%             e = norm_res{2+k};
%             for l = 1:numel(e)
%                 for b = 1:numel(e{l})
%                     nres = nres + norm(e{l}{b})^2;
%                 end
%             end
%         end
%         nres = sqrt(nres);
%         fprintf('|norm_res| = %e \n\n', nres);
%         nv2 = 0;
%         nproj = 0;
%         ntblock = 0;
%         for k = 1:K
%             e = m.v2(k,1);
%             e1 = m.proj(k,1);
%             e2 = m.t_block(k,1);
%             for l = 1:numel(e{1})
%                 for b = 1:numel(e{1}{l})
%                     nv2 = nv2 + norm(e{1}{l}{b}, 'fro');
%                     nproj = nproj + norm(e1{1}{l}{b}, 'fro');
%                     ntblock = ntblock + norm(e2{1}{l}{b}, 'fro');
%                 end
%             end
%         end
%         fprintf('|v2| = %e, |proj| = %e, |t_block| = %e \n\n', nv2, nproj, ntblock);
%         fprintf('|g| = %e \n\n', norm(g(:)))
%         param
        %!--
        %%
    
        % reweight_step_count = reweight_step_count + 1;
        % reweight_last_step_iter = t;
        % rw_counts = rw_counts + 1;   

        if (reweight_step_count >= param.total_reweights)
            % param.reweight_max_reweight_itr = t+1;
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
param.init_v0 = v0_{1};
v0_{1} = [];
param.init_v1 = v1_{2};
v1_{2} = [];
param.init_weights0 = weights0_{1};
weights0_{1} = [];
param.init_weights1 = weights1_{2};
weights1_{2} = [];

% to be completely modified (within spmd function?)
% Calculate residual images
res = zeros(size(xsol));
spmd
    if labindex > 2
        res_ = compute_residual_images(xsol(:,:,c_chunks{labindex-2}), yp, Gp, Ap, Atp, Wp);
    end
end
m = matfile([name, '_', ...
            num2str(param.ind) '_' num2str(param.gamma) '_' num2str(reweight_step_count) '.mat'], ...
            'Writable', true);
m.param = param;
m.res = zeros(size(xsol));
m.g = g;
m.epsilon = cell(K, 1);
m.v2 = cell(K, 1);
m.proj = cell(K, 1);
m.t_block = cell(K, 1);
m.norm_res = cell(K, 1);
m.v0 = v0_{1};
m.v1 = v1_{2};
m.weights0 = weights0_{1};
m.weights1 = weights1_{2};
% Retrieve variables from workers
% data nodes
for k = 1:K
    res(:,:,c_chunks{k}) = res_{2+k};
    m.res(:,:,c_chunks{k}) = res_{2+k};
    res_{2+k} = [];
    m.v2(k,1) = v2_(2+k);
    m.proj(k,1) = proj_(2+k);
    m.t_block(k,1) = t_block(2+k);
    m.epsilon(k,1) = epsilonp(2+k);
    m.norm_res(k,1) = norm_res(2+k);
end
m.xsol = xsol;
epsilon = m.epsilon; % see if necessary
norm_res_out = sqrt(sum(sum(sum((m.res).^2))));

% Update param structure and save
param.reweight_alpha = reweight_alpha;
param.init_reweight_step_count = reweight_step_count;
param.init_reweight_last_iter_step = t;
param.init_t_start = t+1;
m.param = param;

% SNR (computed only on the master node)
sol = reshape(xsol(:),numel(xsol(:))/c,c);
SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
psnrh = zeros(c,1);
for i = 1:c
    psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
end
SNR_average = mean(psnrh);
m.SNR = SNR;
m.SNR_average = SNR_average;
m.end_iter = end_iter;
m.t_master = t_master;
m.t_l21 = t_l21;
m.t_nuclear = t_nuclear;
m.t_data = t_data;
m.rel_val = rel_val;
clear m

% Final log
% nuclear norm
nuclear = nuclear_norm(xsol);
% l21 norm
l21 = compute_sara_prior(xsol, Psit, s);
if (param.verbose > 0)
    if (flag == 1)
        fprintf('Solution found\n');
        fprintf('Iter %i\n',t);
        fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
        fprintf('epsilon = %e, residual = %e\n', norm_epsilon_check,norm_residual_check);
        fprintf('SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf('Iter %i\n',t);
        fprintf('N-norm = %e, L21-norm = %e, rel_val = %e\n', nuclear, l21, rel_val(t));
        fprintf('epsilon = %e, residual = %e\n', norm_epsilon_check,norm_residual_check);
        fprintf('SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
    end
end

end_iter = end_iter(end_iter > 0);
rel_val = rel_val(1:numel(end_iter));

end
