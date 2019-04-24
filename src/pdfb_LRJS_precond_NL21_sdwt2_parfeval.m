function [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = pdfb_LRJS_precond_NL21_sdwt2_parfeval(y, epsilon, A, At, pU, G, W, Sp, Spt, param, X0, Qx, Qy, K, wavelet, L, nlevel, c_chunks, Psit)

%PARFEVAL version: use parfeval for all the priors, deal with the data 
% fidelity term in a single place (one single block)
% cannot recover the time taken by one step with the other one... to be
% investigated further to get a meaningful code).

% This function solves:
%
% min || X ||_* + lambda * ||Psit(X)||_2,1   s.t.  || Y - A(X) ||_2 <= epsilon and x>=0
%
% Author: Abdullah Abdulaziz

% Useful functions for the projection
% sc = @(z,eps) z*min(eps/norm(z(:)), 1); % scaling
hardt = @(z) max(real(z), 0); %thresholding negative values

c = size(y,2);

% [P.-A.] number of wavelet dictionaries
P = length(wavelet);

% oversampling vectorized data length
No = size(W{1}{1}, 1);

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
for qx = 1:Qx % can be created in parallel if needed
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        I(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1];
        dims(q, :) = [rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
    end
end
clear rg_y rg_x;

%%- begin initialization sdwt2
% [P.-A.] instantiate auxiliary variables for sdwt2
[~, ~, I_overlap_ref, dims_overlap_ref, I_overlap, dims_overlap, ...
    ~, ~, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices([M, N], I, dims, nlevel, wavelet, L);
id_dirac = find(ismember(wavelet, 'self'), 1);
dirac_present = ~isempty(id_dirac);

% [P.-A.]
% define parpool for the parallalization and definition of Composite
% objects (see where this can be defined)
numworkers = K*Q; % total number of workers... order: all the K per q (K = row index)
cirrus_cluster = parcluster('local');
cirrus_cluster.NumWorkers = numworkers;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
if cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads > ncores
    exit(1);
end

% [P.-A.] /!\ only simple indexing allowed into Composite objects from the
% master
% instantiate l21 dual variables
parpool(cirrus_cluster, numworkers); % override default preference

% define composite variables (local to a given worker)
% dimension of the ghost cells
overlap_g_south = cell(numworkers, 1);
overlap_g_east = cell(numworkers, 1);
overlap_g_south_east = cell(numworkers, 1);
overlap = cell(numworkers, 1);
% initialize composite variables and constants
% amount of overlap of the neighbour (necessary to define the ghost cells properly)
for q = 1:Q
    % additional composite variables (for the zero padding, see if fewer elements can be used)
    overlap{q} = max(dims_overlap{q}) - dims(q,:); % amount of overlap necessary for each facet
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

%Initial dual variables
%% to be completely updated (add options for warm-restart later on)
% [P.-A.] l21 dual variables
v0 = cell(numworkers, 1);
weights0 = cell(numworkers, 1);
v1 = cell(numworkers, 1);
weights1 = cell(numworkers, 1);
for q = 1:Q
    p = prod(Ncoefs{q}, 2);
    if dirac_present
        sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) + prod(dims(q,:));
    else
        sz = 3*sum(p) - 2*sum(p(nlevel+1:nlevel+1:end));
    end
    
    for i = 1:K
        qi = (q-1)*K + i;
        % nuclear-norm prior
        v0{qi} = zeros(prod(dims_overlap_ref(q,:)), c_chunks(i));
        weights0{qi} = ones(min(prod(dims_overlap_ref(q, :)), c_chunks(i)), 1);
        % joint-sparsity prior
        v1{qi} = zeros(sz, c_chunks(i));
        weights1{qi} = ones(sz, c_chunks(i));
    end
end
clear sz
%%

if isfield(param,'init_v2')
    v2 = param.init_v2;
    for i = 1 : c
        u2{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            u2{i}{j} = zeros(size(G{i}{j}, 2), 1);
        end
    end
    r2 = v2;
    fprintf('v2 uploaded \n\n')
else
    for i = 1 : c
        v2{i} = cell(length(G{i}),1);
        u2{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            v2{i}{j} = zeros(length(y{i}{j}) ,1);
            u2{i}{j} = zeros(size(G{i}{j}, 2), 1);
        end
    end
    r2 = v2;
    fprintf('v2 NOT uploaded \n\n')
end


%Initial primal gradient
if isfield(param,'init_g')
    g = param.init_g;
    fprintf('g uploaded \n\n')
else
    g = zeros(size(xsol));
    fprintf('g NOT uploaded \n\n')
end

Ftx = zeros(size(xsol));


% Initialise projection
if isfield(param,'init_proj')
    proj = param.init_proj;
    fprintf('proj uploaded \n\n')
else
    proj = cell(c, 1);
    for i = 1 : c
        proj{i} = cell(length(G{i}),1);
        Fx = A(xsol(:,:,i));
        for j = 1 : length(G{i})
            r2{i}{j} = G{i}{j} * Fx(W{i}{j});
            [proj{i}{j}, ~] = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, r2{i}{j}, y{i}{j}, pU{i}{j}, epsilon{i}{j}, zeros(size(y{i}{j})), param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);
        end
    end
    fprintf('proj NOT uploaded \n\n')
end


if isfield(param,'init_t_block')
    t_block = param.init_t_block;
    t_start = param.init_t+1;
    reweight_last_step_iter = param.init_t;
    reweight_step_count = param.reweight_step_count+1;
    % rw_counts is an index for the reweight steps vector
    rw_counts = 1;
    fprintf('t t_block uploaded \n\n')
else
    t_block = cell(c, 1);
    for i = 1 : c
        t_block{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            t_block{i}{j} = 0;
        end
    end
    t_start = 1;
    reweight_last_step_iter = 0;
    reweight_step_count = 0;
    rw_counts = 1;
    fprintf('t t_block NOT uploaded \n\n')
end

count_eps_update_down = 0;
count_eps_update_up = 0;

reweight_alpha = param.reweight_alpha;
reweight_alpha_ff = param.reweight_alpha_ff;
reweight_steps = param.reweight_steps;

%Step sizes computation

%Step size for the dual variables
sigma0 = 1.0/param.nu0;
sigma1 = 1.0/param.nu1;
sigma2 = 1.0/param.nu2;

%Step size primal
tau = 0.99/(sigma0*param.nu0 + sigma1*param.nu1 + sigma2*param.nu2);

sigma00 = tau*sigma0;
sigma11 = tau*sigma1;
sigma22 = tau*sigma2;

flag = 0;

beta0 = param.gamma0/sigma0;
beta1 = param.gamma/sigma1;

Gt = cell(size(G));
for i = 1 : c
    Gt{i} = cell(length(G{i}),1);
    for j = 1 : length(G{i})
        Gt{i}{j} = G{i}{j}';
    end
end

A = afclean(A);
At = afclean(At);


% Main loop. Sequential.
% util_create_pool(param.num_workers);
% maxNumCompThreads(param.num_workers);

start_loop = tic;

for t = t_start : param.max_iter
    
    %fprintf('Iter %i\n',t);
    start_iter = tic;
    
    %% Primal update
    prev_xsol = xsol;
    xsol = hardt(xsol - g);
    xhat = 2*xsol - prev_xsol;
    xhat_split = Sp(xhat);
    
    %% Relative change of objective function
    rel_fval(t) = norm(xsol(:) - prev_xsol(:))/norm(xsol(:));
    % Free memory
    %prev_xsol = [];
    
    %% Dual variables update (spmd)
    
    %% Nuclear norm / L-2,1 function update
    
    % communicate x_overlap (full communications) [sub-optimal here]
    for qi = 1:numworkers
        [i, q] = ind2sub([K, Q], qi);
        x_overlap = xhat_split{i}(I_overlap_ref(q, 1)+1:I_overlap_ref(q, 1)+dims_overlap_ref(q, 1), ...
            I_overlap_ref(q, 2)+1:I_overlap_ref(q, 2)+dims_overlap_ref(q, 2), :);
        
        s(qi) = parfeval(@update_prior, 3, x_overlap, v0{qi}, v1{qi}, ...
                weights0{qi}, weights1{qi}, beta0, beta1, ...
                dims_overlap_ref(q, :), offset, status(q, :), offsetL(q,:), offsetR(q,:), ...
                I(q, :), dims(q,:), I_overlap{q}, dims_overlap{q}, Ncoefs{q}, ...
                nlevel, wavelet, temLIdxs{q}, temRIdxs{q}, sigma00, sigma11); % to be debugged
    end
    
    %% L2 ball projection update (currently done on a single node: is there any point in having multiple data blocks then?)
    % this step takes much more time than the step above...
    counter = 1;
    for i = 1 : c
        Fx = A(xhat(:,:,i));
        g2 = zeros(No,1);
        for j = 1 : length(G{i})
            r2{i}{j} = G{i}{j} * Fx(W{i}{j});
            [proj{i}{j}, ~] = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, r2{i}{j}, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);
            v2{i}{j} = v2{i}{j} + pU{i}{j} .* r2{i}{j} - pU{i}{j} .* proj{i}{j};
            
            u2{i}{j} = Gt{i}{j} * v2{i}{j};
            g2(W{i}{j}) = g2(W{i}{j}) + u2{i}{j};
            
            % norm of residual
            norm_res{i}{j} = norm(r2{i}{j} - y{i}{j});
            residual_check(counter) = norm_res{i}{j};
            epsilon_check(counter) = epsilon{i}{j};
            counter = counter + 1;
        end
        Ftx(:,:,i) = real(At(g2));
    end
    % Free memory
    %g2=[]; Fx=[];
    
    %% Update primal gradient (part 1)
    go = cell(numworkers, 1);
    for qi = 1 : numworkers
        [idx, v0_, v1_, go_] = fetchNext(s);
        v0{idx} = v0_;
        v1{idx} = v1_;
        go{idx} = go_;
    end
    
    g_ = cell(K, 1);
    for i = 1:K
        g_{i} = zeros(M, N, c_chunks(i));
        for q = 1:Q
            qi = (q-1)*K + i;
            % aggregate elements back in place (overlap)
            g_{i} = place2DSegment(g_{i}, go{qi}, I_overlap_ref(q, :), ...
                    dims_overlap_ref(q, :));
        end
    end
    
    %% Update primal gradient (part 2)
    
    g = Spt(g_) + sigma22*Ftx;
    % Free memory
    %g0=[]; g1=[]; Ftx=[];
    
    end_iter(t) = toc(start_iter);
    fprintf('Iter = %i, Time = %e\n',t,end_iter(t));
    
    %% Display
    if ~mod(t,100)
        
        % [P.-A.] to be modified (facets are missing...)
        xhatm = reshape(xsol,numel(xsol)/c,c);
        [~,S0,~] = svd(xhatm,'econ');
        nuclear = norm(diag(S0),1); % doesn't make sense to compute this here
        
        % [P.-A.] to be modified (idem, facets are missing)
        l2 = sqrt(sum(abs(Psit(xsol)).^2,2));
        l21 = norm(l2(:),1);
        
        %SNR
        sol = reshape(xsol(:),numel(xsol(:))/c,c);
        SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
        psnrh = zeros(c,1);
        for i = 1:c
            psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
        end
        SNR_average = mean(psnrh);
        
        %Log
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
            fprintf(' epsilon = %e, residual = %e\n', norm(epsilon_check),norm(residual_check));
            fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
        end
        
        figure(1),
        subplot(2,2,1);
        imagesc(log10(max(flip(x0(:,:,1)),0))); hold on; colorbar; axis image; axis off; colormap(cubehelix); caxis([-3.5, 0]);
        subplot(2,2,2);
        imagesc(log10(max(flip(xsol(:,:,1)),0))); hold on; colorbar; axis image; axis off; colormap(cubehelix); caxis([-3.5, 0]);
        subplot(2,2,3);
        imagesc(log10(max(flip(x0(:,:,end)),0))); hold on; colorbar; axis image; axis off; colormap(cubehelix); caxis([-3.5, 0]);
        subplot(2,2,4);
        imagesc(log10(max(flip(xsol(:,:,end)),0))); hold on; colorbar; axis image; axis off; colormap(cubehelix); caxis([-3.5, 0]);
        pause(0.1)
        
        fitswrite(x0,'./x0.fits');
        fitswrite(xsol,'./xsol.fits');
    end
    
    %% Save
    if ~mod(t,200000)
        
        % Calculate residual images:
        for i = 1 : c
            Fx = A(xsol(:,:,i));
            g2 = zeros(No,1);
            for j = 1 : length(G{i})
                res_f{i}{j} = y{i}{j} - G{i}{j} * Fx(W{i}{j});
                u2{i}{j} = Gt{i}{j} * res_f{i}{j};
                g2(W{i}{j}) = g2(W{i}{j}) + u2{i}{j};
            end
            res(:,:,i) = real(At(g2));
        end
        
        fitswrite(xsol,['./reweight_res/xsol_5G_',num2str(t),'.fits']);
        fitswrite(res,['./reweight_res/res_5G_',num2str(t),'.fits']);
        
    end
    
    %% Global stopping criteria
    %     prod(prod(residual_check < param.adapt_eps_tol_out*epsilon_check)) && prod(prod(residual_check > param.adapt_eps_tol_in*epsilon_check))
    if t>1 && rel_fval(t) < param.rel_obj && reweight_step_count > param.total_reweights && ...
            (norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check))
        flag = 1;
        break;
    end
    
    %% Update epsilons
    if param.use_adapt_eps && t > param.adapt_eps_start
        for i = 1 : c
            for  j = 1 : length(G{i})
                if  norm_res{i}{j} < param.adapt_eps_tol_in * epsilon{i}{j}
                    if t > t_block{i}{j} + param.adapt_eps_steps && rel_fval(t) < param.adapt_eps_rel_obj
                        epsilon{i}{j} = norm_res{i}{j} + (-norm_res{i}{j} + epsilon{i}{j}) * (1 - param.adapt_eps_change_percentage);
                        t_block{i}{j} = t;
                        count_eps_update_down = count_eps_update_down + 1;
                        fprintf('Updated  epsilon DOWN: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j},norm_res{i}{j},j,i);
                    end
                end
                
                if  norm_res{i}{j} > param.adapt_eps_tol_out * epsilon{i}{j}
                    if t > t_block{i}{j} + param.adapt_eps_steps && rel_fval(t) < param.adapt_eps_rel_obj
                        epsilon{i}{j} = epsilon{i}{j} + (norm_res{i}{j} - epsilon{i}{j}) * param.adapt_eps_change_percentage;
                        t_block{i}{j} = t;
                        count_eps_update_up = count_eps_update_up + 1;
                        fprintf('Updated  epsilon UP: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j},norm_res{i}{j},j,i);
                    end
                end
            end
        end
    end
    
    %% Reweighting
    if (param.step_flag && rel_fval(t) < param.reweight_rel_obj) % && (norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check))
        reweight_steps = [t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size)];
        param.step_flag = 0;
    end
    
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
            (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
            t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);
        
        xsol_split = Sp(xsol);
        
        for qi = 1:numworkers % ti be done in a parfeval / parfor or spmd
            [i, q] = ind2sub([K, Q], qi);
            x_sol_overlap = xsol_split{i}(I_overlap_ref(q, 1)+1:I_overlap_ref(q, 1)+dims_overlap_ref(q, 1), ...
                I_overlap_ref(q, 2)+1:I_overlap_ref(q, 2)+dims_overlap_ref(q, 2), :);
            
            % reweight nuclear norm
            sol = reshape(x_sol_overlap, [size(x_sol_overlap, 1)*size(x_sol_overlap, 2), c_chunk(i)]);
            [~,S00,~] = svd(sol,'econ');
            d_val0 = abs(diag(S00));
            weights0{qi} = reweight_alpha ./ (reweight_alpha + d_val0);
            weights0{qi}(d_val0 > max(d_val0) * param.reweight_abs_of_max) = 0;
            
            % reweight l21 norm
            zerosNum = dims_overlap_ref(q, :) + offsetL(q, :) + offsetR(q, :); % offset for the zero-padding (to be checked again...)
            x_ = zeros([zerosNum, size(x_sol_overlap, 3)]);
            x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_sol_overlap;
            w = zeros(size(v1{qi}));
            for l = 1 : size(x_, 3)
                w(:, l) = sdwt2_sara(x_(:, :, l), Iq, dims_q, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q);
            end
            d_val1 = sqrt(sum(abs((w)).^2,2));
            weights1{qi} = reweight_alpha ./ (reweight_alpha + d_val1);
        end
        reweight_alpha = reweight_alpha_ff .* reweight_alpha;
        
        %for i = 1 : K
        %    sol = reshape(xsol_split{i},numel(xsol_split{i})/c_chunk(i),c_chunk(i));
        %    [~,S00,~] = svd(sol,'econ');
        %    d_val0 = abs(diag(S00));
        %    weights0{i} = reweight_alpha ./ (reweight_alpha + d_val0);
        %    weights0{i}(d_val0 > max(d_val0) * param.reweight_abs_of_max) = 0;
        %end
        
        % [P.-A.]
        %for i = 1 : K
        %    d_val1 = sqrt(sum(abs(Psit(xsol_split{i})).^2,2));
        %    weights1{i} = reweight_alpha ./ (reweight_alpha + d_val1);
        %    weights1{i}(d_val1 > max(d_val1) * param.reweight_abs_of_max) = 0;
        %end
        %reweight_alpha = reweight_alpha_ff .* reweight_alpha;
        %-- end modifications
        
        if (reweight_step_count >= param.total_reweights)
            param.reweight_max_reweight_itr = t+1;
            fprintf('\n\n No more reweights \n\n');
            break;
        end
        
        % Calculate residual images:
        for i = 1 : c
            Fx = A(xsol(:,:,i));
            g2 = zeros(No,1);
            for j = 1 : length(G{i})
                res_f{i}{j} = y{i}{j} - G{i}{j} * Fx(W{i}{j});
                u2{i}{j} = Gt{i}{j} * res_f{i}{j};
                g2(W{i}{j}) = g2(W{i}{j}) + u2{i}{j};
            end
            res(:,:,i) = real(At(g2));
        end
        
        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;
        rw_counts = rw_counts + 1;
        
    end
    
end
end_loop = toc(start_loop)

% Calculate residual images:
for i = 1 : c
    Fx = A(xsol(:,:,i));
    g2 = zeros(No,1);
    for j = 1 : length(G{i})
        res_f{i}{j} = y{i}{j} - G{i}{j} * Fx(W{i}{j});
        u2{i}{j} = Gt{i}{j} * res_f{i}{j};
        g2(W{i}{j}) = g2(W{i}{j}) + u2{i}{j};
    end
    res(:,:,i) = real(At(g2));
end

%Final log

xhatm = reshape(xsol,numel(xsol)/c,c);
[~,S0,~] = svd(xhatm,'econ');
nuclear = norm(diag(S0),1);

l2 = sqrt(sum(abs(Psit(xsol)).^2,2));
l21 = norm(l2(:),1);

%SNR
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
        fprintf(' epsilon = %e, residual = %e\n', norm(epsilon_check),norm(residual_check));
        fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf('Iter %i\n',t);
        fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
        fprintf(' epsilon = %e, residual = %e\n', norm(epsilon_check),norm(residual_check));
        fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
    end
end

end_iter = end_iter(end_iter > 0);
rel_fval = rel_fval(1:numel(end_iter));

end

% function [v0, g0] = run_par_nuclear(v0, xhat, weights0, beta0)
% 
% [M,N,c] = size(xhat);
% xhatm = reshape(xhat,numel(xhat)/c,c);
% [U0,S0,V0] = svd(v0 + xhatm,'econ');
% v0 = v0 + xhatm - (U0*diag(max(diag(S0) - beta0 * weights0, 0))*V0');
% g0 = reshape(v0,M,N,c);
% %nuclear(t) = norm(diag(S0),1);
% 
% end

% function [v1, u1] = run_par_l21(v1, Psit, Psi, xhat, weights1, beta1)
%
%     r1 = v1 + Psit(xhat);
%     l2 = sqrt(sum(abs(r1).^2,2));
%     l2_soft = max(l2 - beta1*weights1, 0)./ (l2+eps);
%     v1 = r1 - (l2_soft .* r1);
%     u1 = Psi(v1);
%     %l21 = norm(l2(:),1);
% end

function [v0, v1, g] = update_prior(xhat, v0, v1, weights0, weights1, ...
    beta0, beta1, dims_overlap_ref_q, offset, status_q, offsetLq, offsetRq, ...
    Iq, dims_q, I_overlap_q, dims_overlap_q, Ncoefs_q, nlevel, wavelet, ...
    temLIdxs_q, temRIdxs_q, sigma00, sigma11)
% see if I can avoid communicating all these informations each time (quite cumbersome...)

% nuclear norm prior (see if I remove the overlap here...)
[M,N,c] = size(xhat);
xhatm = reshape(xhat,M*N,c);
[U0,S0,V0] = svd(v0 + xhatm,'econ');
v0 = v0 + xhatm - (U0*diag(max(diag(S0) - beta0 * weights0, 0))*V0');
% g0 = reshape(v0,M,N,c);
%nuclear(t) = norm(diag(S0),1);

% l21 prior
zerosNum = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding (to be checked again...)
x_ = zeros([zerosNum, size(xhat, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = xhat;
g1 = zeros(size(xhat));

l21_ = 0;
for l = 1 : size(x_, 3)
    u_old = v1(:, l);
    w = sdwt2_sara(x_(:, :, l), Iq, dims_q, offset, status_q, nlevel, wavelet, Ncoefs_q);
    w = v1(:, l) +  w;
    l2 = sqrt(sum(abs(w).^2,2));
    l2_soft = max(l2 - beta1*weights1(:, l), 0)./(l2+eps);
    v1(:, l) = w - l2_soft.*w;
    
    g1(:, :, l) = isdwt2_sara(v1(:, l)-u_old, Iq, dims_q, I_overlap_q, dims_overlap_q, Ncoefs_q, nlevel, wavelet, temLIdxs_q, temRIdxs_q);
    
    % local L21 norm of current solution
    l21_ = l21_ + norm(l2(:),1);
end

g = sigma00*reshape(v0,M,N,c) + sigma11*g1;

end
