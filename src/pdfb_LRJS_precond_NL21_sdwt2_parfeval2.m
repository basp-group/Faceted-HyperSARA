function [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = pdfb_LRJS_precond_NL21_sdwt2_parfeval2(y, epsilon, A, At, pU, G, W, param, X0, Qx, Qy, wavelet, L, nlevel, Psit)

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

% oversampling vectorized data length
No = size(W{1}{1}, 1);

% number of pixels
[M, N] = size(At(zeros(No, 1)));

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
[~, ~, I_overlap_ref, dims_overlap_ref, I_overlap, dims_overlap, ...
    ~, ~, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices([M, N], I, dims, nlevel, wavelet, L);
id_dirac = find(ismember(wavelet, 'self'), 1);
dirac_present = ~isempty(id_dirac);

% [P.-A.]
% define parpool for the parallalization and definition of Composite
% objects (see where this can be defined)
numworkers = 2*Q; % total number of workers (faceted nuclear and l21 norms)
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
overlap_g_south = cell(Q, 1);
overlap_g_east = cell(Q, 1);
overlap_g_south_east = cell(Q, 1);
overlap = cell(Q, 1);
% initialize composite variables and constants
% amount of overlap of the neighbour (necessary to define the ghost cells properly)
for q = 1:Q
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
v0 = cell(Q, 1);
weights0 = cell(Q, 1);
v1 = cell(Q, 1);
weights1 = cell(Q, 1);
for q = 1:Q
    p = prod(Ncoefs{q}, 2);
    if dirac_present
        sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) + prod(dims(q,:));
    else
        sz = 3*sum(p) - 2*sum(p(nlevel+1:nlevel+1:end));
    end
    
    % nuclear-norm prior
    v0{q} = zeros(prod(dims(q,:)), c);
    weights0{q} = zeros(min(prod(dims(q, :)), c), 1);
    % joint-sparsity prior
    v1{q} = zeros(sz, c);
    weights1{q} = zeros(sz, c);
end
clear sz
%%

if isfield(param,'init_v2')
    v2 = param.init_v2;
    r2 = v2;
    fprintf('v2 uploaded \n\n')
else
    v2 = cell(c, 1);
    for i = 1 : c
        v2{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            v2{i}{j} = zeros(length(y{i}{j}) ,1);
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
    
    %% Relative change of objective function
    rel_fval(t) = norm(xsol(:) - prev_xsol(:))/norm(xsol(:));
    % Free memory
    %prev_xsol = [];
    
    %% Dual variables update (parfeval)
    
    %% Nuclear norm function update
    for q = 1:Q
         s(q) = parfeval(@run_par_nuclear, 2, v0{q}, xhat(I(q, 1)+1:I(q, 1)+dims(q, 1), ...
                I(q, 2)+1:I(q, 2)+dims(q, 2), :), weights0{q}, beta0);
    end  
    
    %% L-2,1 function update
    for q = 1:Q
        x_overlap = xhat(I_overlap_ref(q, 1)+1:I_overlap_ref(q, 1)+dims_overlap_ref(q, 1), ...
            I_overlap_ref(q, 2)+1:I_overlap_ref(q, 2)+dims_overlap_ref(q, 2), :);
        zerosNum = dims_overlap_ref(q, :) + offsetL(q, :) + offsetR(q, :);
        
        f(q) = parfeval(@run_par_l21, 2, v1{q}, I(q, :), dims(q, :), I_overlap{q}, ...
            dims_overlap{q}, offset, status(q, :), nlevel, wavelet, Ncoefs{q}, ...
            zerosNum, temLIdxs{q}, temRIdxs{q}, offsetL(q, :), offsetR(q, :), x_overlap, weights1{q}, beta1);
    end         
    
    %% L2 ball projection update
    counter = 1;
    for i = 1 : c
        Fx = A(xhat(:,:,i));
        g2 = zeros(No,1);
        for j = 1 : length(G{i})
            r2{i}{j} = G{i}{j} * Fx(W{i}{j});
            [proj{i}{j}, ~] = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, r2{i}{j}, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);
            v2{i}{j} = v2{i}{j} + pU{i}{j} .* r2{i}{j} - pU{i}{j} .* proj{i}{j};
            
            u2 = Gt{i}{j} * v2{i}{j};
            g2(W{i}{j}) = g2(W{i}{j}) + u2;
            
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
    g = zeros(size(xsol)); % no spectral faceting
    v0 = cell(Q, 1);
    for q = 1 : Q
        [idx, v0_, g0_] = fetchNext(s);
        v0{idx} = v0_;
        g(I(idx, 1)+1:I(idx, 1)+dims(idx, 1), I(idx, 2)+1:I(idx, 2)+dims(idx, 2), :) = ...
        g(I(idx, 1)+1:I(idx, 1)+dims(idx, 1), I(idx, 2)+1:I(idx, 2)+dims(idx, 2), :) + g0_;
    end
    g = sigma00*g;
    
    for q = 1 : Q
        [idx, v1_, g1_] = fetchNext(f);
        v1{idx} = v1_;
        g = place2DSegment(g, sigma11*g1_, I_overlap_ref(idx, :), ...
                    dims_overlap_ref(idx, :));
    end
    g = g + sigma22*Ftx;
    % Free memory
    %g0=[]; g1=[]; Ftx=[];
    
    end_iter(t) = toc(start_iter);
    fprintf('Iter = %i, Time = %e\n',t,end_iter(t));
    
    %% Display
    if ~mod(t,100)
        
        % [P.-A.] to be modified (facets are missing...)
        xhatm = reshape(xsol,numel(xsol)/c,c);
        [~,S0,~] = svd(xhatm,'econ');
        nuclear = norm(diag(S0),1);
        
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

%         fitswrite(x0,'./x0.fits');
%         fitswrite(xsol,'./xsol.fits');
    end
    
    %% Save
    if ~mod(t,200000)
        
        % Calculate residual images:
        for i = 1 : c
            Fx = A(xsol(:,:,i));
            g2 = zeros(No,1);
            for j = 1 : length(G{i})
                res_f = y{i}{j} - G{i}{j} * Fx(W{i}{j});
                u2 = Gt{i}{j} * res_f{i}{j};
                g2(W{i}{j}) = g2(W{i}{j}) + u2;
            end
            res(:,:,i) = real(At(g2));
        end
        
        fitswrite(xsol,['./reweight_res/xsol_5G_',num2str(t),'.fits']);
        fitswrite(res,['./reweight_res/res_5G_',num2str(t),'.fits']);
        
    end
    
    %% Global stopping criteria
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
    if (param.step_flag && rel_fval(t) < param.reweight_rel_obj)
        reweight_steps = [t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size)];
        param.step_flag = 0;
    end
    
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
            (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
            t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);
        
        for q = 1:Q            
            sol = reshape(xsol(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :),prod(dims(q, :)),c);
            [~,S00,~] = svd(sol,'econ');
            d_val0 = abs(diag(S00));
            weights0{q} = reweight_alpha ./ (reweight_alpha + d_val0);
            weights0{q}(d_val0 > max(d_val0) * param.reweight_abs_of_max) = 0;
            
            zerosNum = dims_overlap_ref(q, :) + offsetL(q, :) + offsetR(q, :); % offset for the zero-padding (to be checked again...)
            x_ = zeros([zerosNum, c]);
            x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = ...
                xsol(I_overlap_ref(q, 1)+1:I_overlap_ref(q, 1)+dims_overlap_ref(q, 1), I_overlap_ref(q, 2)+1:I_overlap_ref(q, 2)+dims_overlap_ref(q, 2), :);
            w = zeros(size(v1{q}));
            for l = 1 : size(x_, 3)
                w(:, l) = sdwt2_sara(x_(:, :, l), I(:, q), dims(q, :), offset, status(q, :), nlevel, wavelet, Ncoefs{q});
            end
            d_val1 = sqrt(sum(abs((w)).^2,2));
            weights1{q} = reweight_alpha ./ (reweight_alpha + d_val1);
            weights1{q}(d_val1 > max(d_val1) * param.reweight_abs_of_max) = 0;
        end
        reweight_alpha = reweight_alpha_ff .* reweight_alpha;
        
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
        res_f = y{i}{j} - G{i}{j} * Fx(W{i}{j});
        u2 = Gt{i}{j} * res_f;
        g2(W{i}{j}) = g2(W{i}{j}) + u2;
    end
    res(:,:,i) = real(At(g2));
end

% Final log
xhatm = reshape(xsol,numel(xsol)/c,c);
[~,S0,~] = svd(xhatm,'econ');
nuclear = norm(diag(S0),1);

l2 = sqrt(sum(abs(Psit(xsol)).^2,2));
l21 = norm(l2(:),1);

% SNR
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


end

function [v0, g0] = run_par_nuclear(v0, xhat, weights0, beta0)

[M,N,c] = size(xhat);
xhatm = reshape(xhat,numel(xhat)/c,c);
[U0,S0,V0] = svd(v0 + xhatm,'econ');
v0 = v0 + xhatm - (U0*diag(max(diag(S0) - beta0 * weights0, 0))*V0');
g0 = reshape(v0,M,N,c);

end

function [v1, g1] = run_par_l21(v1, Iq, dims_q, I_overlap_q, dims_overlap_q, offset, status_q, nlevel, ...
    wavelet, Ncoefs_q, zerosNum, temLIdxs_q, temRIdxs_q, offsetLq, offsetRq, xhat, weights1, beta1)

%zerosNum = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding (to be checked again...)
x_ = zeros([zerosNum, size(xhat, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = xhat;
g1 = zeros(size(xhat));

for l = 1 : size(x_, 3)
    w = sdwt2_sara(x_(:, :, l), Iq, dims_q, offset, status_q, nlevel, wavelet, Ncoefs_q);
    w = v1(:, l) +  w;
    l2 = sqrt(sum(abs(w).^2,2));
    l2_soft = max(l2 - beta1*weights1(:, l), 0)./(l2+eps);
    v1(:, l) = w - l2_soft.*w;
    
    g1(:, :, l) = isdwt2_sara(v1(:, l), Iq, dims_q, I_overlap_q, dims_overlap_q, Ncoefs_q, nlevel, wavelet, temLIdxs_q, temRIdxs_q);
end

end
