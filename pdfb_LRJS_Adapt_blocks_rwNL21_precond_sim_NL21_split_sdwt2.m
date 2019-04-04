function [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res] = pdfb_LRJS_Adapt_blocks_rwNL21_precond_sim_NL21_split_sdwt2(y, epsilon, A, At, pU, G, W, Sp, Spt, param, X0, Qx, Qy, Qc, wavelet, L, nlevel, c_chunks, spectral_facets)

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

% useless
% for i = 1 : c
%     x0(:,:,i) = reshape(X0(:,i),M,N);
% end
x0 = reshape(X0, [M, N, c]);

% [P.-A.]
% define spatial facets (no overlap)
Q = Qx*Qy;
rg_y = domain_decomposition(Qy, N(1));
rg_x = domain_decomposition(Qx, N(2));
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
    ~, ~, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices(N, I, dims, nlevel, wavelet, L);
id_dirac = find(ismember(wavelet, 'self'), 1);
dirac_present = ~isempty(id_dirac);

% [P.-A.]
% define parpool for the parallalization and definition of Composite
% objects (see where this can be defined)
numworkers = Qc*Q; % total number of workers... order: all the Qc per q (Qc = row index)
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

% define parallel constants(known by each worker)
Qyp = parallel.pool.Constant(Qy);
Qxp = parallel.pool.Constant(Qx);
Np = parallel.pool.Constant([M, N]);
waveletp = parallel.pool.Constant(wavelet);
nlevelp = parallel.pool.Constant(nlevel);
offsetp = parallel.pool.Constant(offset);

% define composite variables (local to a given worker)
x_overlap = Composite(numworkers);
x_sol_overlap = Composite(numworkers);
% wavelet auxiliary variables
Iq = Composite(numworkers);
dims_q = Composite(numworkers);
temLIdxs_q = Composite(numworkers);
temRIdxs_q = Composite(numworkers);
I_overlap_q = Composite(numworkers);
dims_overlap_q = Composite(numworkers);
dims_overlap_ref_q = Composite(numworkers);
status_q = Composite(numworkers);
Ncoefs_q = Composite(numworkers);
offsetLq = Composite(numworkers);
offsetRq = Composite(numworkers);
% dimension of the ghost cells
overlap_g_south = Composite(numworkers);
overlap_g_east = Composite(numworkers);
overlap_g_south_east = Composite(numworkers);
overlap = Composite(numworkers);
% initialize composite variables and constants
for q = 1:Q
    for i = 1:Qc
        qi = (q-1)*Q + i;
        Iq{qi} = I(q, :);
        dims_q{qi} = dims(q, :);
        temLIdxs_q{qi} = temLIdxs{q};
        temRIdxs_q{qi} = temRIdxs{q};
        I_overlap_q{qi} = I_overlap{q};
        dims_overlap_q{qi} = dims_overlap{q};
        status_q{qi} = status(q, :);
        Ncoefs_q{qi} = Ncoefs{q};
        
        % additional composite variables (for the zero padding, see if fewer elements can be used)
        dims_overlap_ref_q{qi} = dims_overlap_ref(q,:);
        offsetLq{qi} = offsetL(q,:);
        offsetRq{qi} = offsetR(q,:);
        overlap{qi} = max(dims_overlap{q}) - dims(q,:); % amount of overlap necessary for each facet        
    end
end

% amount of overlap of the neighbour (necessary to define the ghost cells properly)
for q = 1:Q
    [qy, qx] = ind2sub([Qy, Qx], q);
    if qy < Qy
        % S (qy+1, qx)
        for i = 1:Qc
            qi = (q-1)*Qc + i;
            overlap_g_south{qi} = overlap{(qx-1)*Qy + qy+1};
        end
        if qx < Qx
            % SE (qy+1, qx+1)
            for i = 1:Qc
                qi = (q-1)*Qc + i;
                overlap_g_south_east{qi} = overlap{qx*Qy + qy+1};
            end
        else
            for i = 1:Qc
                qi = (q-1)*Qc + i;
                overlap_g_south_east{qi} = [0, 0];
            end
        end
    else
        for i = 1:Qc
            qi = (q-1)*Qc + i;
            overlap_g_south{qi} = [0, 0];
            overlap_g_south_east{qi} = [0, 0];
        end
    end
    if qx < Qx
        % E (qy, qx+1)
        for i = 1:Qc
            qi = (q-1)*Qc + i;
            overlap_g_east{qi} = overlap{qx*Qy + qy};
        end
    else
        for i = 1:Qc
            qi = (q-1)*Qc + i;
            overlap_g_east{qi} = [0, 0];
        end
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
if isfield(param,'init_v0')
    v0 = param.init_v0;
    fprintf('v0 uploaded \n\n')
else
    v0 = cell(Qc, 1);
    for i = 1:Qc
        v0{i} = zeros(M*N, c_chunks(i));
    end
    fprintf('v0 NOT uploaded \n\n')
end

if isfield(param,'init_weights0')
    weights0 = param.init_weights0;
    fprintf('weights0 uploaded \n\n')
else
    weights0 = cell(Qc, 1);
    for i = 1 : Qc
        weights0{i} = ones(c_chunks(i),1);
    end
    fprintf('weights0 NOT uploaded \n\n')
end

%Initial dual variables

%% to be completely updated (add options for warm-restart later on)
% [P.-A.] l21 dual variables 
v1 = Composite(numworkers);
weights1 = Composite(numworkers);
for q = 1:Q
    p = prod(Ncoefs{q}, 2);
    if dirac_present
        s = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) + prod(dims(q,:));
    else
        s = 3*sum(p) - 2*sum(p(nlevel+1:nlevel+1:end));
    end
    
    for i = 1:Qc
        qi = (q-1)*Qc + i;
        v1{qi} = zeros(s, c_chunks(i)); % multi-band case
        weights1{qi} = zeros(s, c_chunks(i));
    end
end

%%

if isfield(param,'init_v2')
    v2 = param.init_v2;
    u2 = cell(c, 1);
    for i = 1 : c
        u2{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            u2{i}{j} = zeros(size(G{i}{j}, 2), 1);
        end
    end
    r2 = v2;
    fprintf('v2 uploaded \n\n')
else
    u2 = cell(c, 1);
    v2 = cell(c, 1);
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

g0 = cell(Qc, 1);
Fx = zeros(No,c);
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
reweight_alphap = Composite(numworkers);
for qi = 1:numworkers
    reweight_alphap{qi} = reweight_alpha;
end
reweight_alpha_ff = parallel.pool.Constant(param.reweight_alpha_ff);
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
beta1 = parallel.pool.Constant(param.gamma/sigma1); % on the workers


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
%maxNumCompThreads(12);
% util_create_pool(24);

for t = t_start : param.max_iter
    
    %fprintf('Iter %i\n',t);
    %tic;
    
    %% Primal update
    prev_xsol = xsol;
    xsol = hardt(xsol - g);
    xhat = 2*xsol - prev_xsol;
    xhat_split = Sp(xhat);
    
    % communicate x_overlap (full communications) [sub-optimal here]
    for qi = 1:numworkers
        [i, q] = ind2sub([Qc, Q], qi);
        x_overlap{qi} = xhat_split(I_overlap_ref(q, 1)+1:I_overlap_ref(q, 1)+dims_overlap_ref(q, 1), ...
        I_overlap_ref(q, 2)+1:I_overlap_ref(q, 2)+dims_overlap_ref(q, 2), spectral_facets{i});
    end
    
    %% Relative change of objective function
    rel_fval(t) = norm(xsol(:) - prev_xsol(:))/norm(xsol(:));
    % Free memory
    prev_xsol = [];
    
    %% Dual variables update
    
    %% Nuclear norm function update % [P.-A.] see how to make this parallelization more efficient
    % parfor
    nuclear(t) = 0;
    for i = 1 : Qc
        temp = xhat_split{i};
        xhatm = reshape(temp(:),numel(temp(:))/c_chunk(i),c_chunk(i));
        [U0,S0,V0] = svd(v0{i} + xhatm,'econ');
        nuclear(t) = nuclear(t) + norm(diag(S0),1);
        v0{i} = v0{i} + xhatm - (U0*diag(max(diag(S0) - beta0 * weights0{i}, 0))*V0');
    end
    % Free memory
    U0=[]; S0=[]; V0=[]; xhatm = [];
    
    
    %% L-2,1 function update [P.-A.] to be completely changed here...
    % spmd (update primal variable here?)
    % for k = 1:P
    %     f(k) = parfeval(@run_par_waverec, 3, v1{k}, Psit{k}, Psi{k}, xhat_split, weights1{k}, beta1,c_chunk);
    % end
    
    % spmd update here, check which constants are needed
    spmd
        [v1, u1, l21_] = run_par_waverec(v1, x_overlap, weights1, beta1, Iq, dims_q, I_overlap_q, dims_overlap_q, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q, temLIdxs_q, temRIdxs_q, offsetLq, offsetRq, dims_overlap_ref_q);
        % reduction with optimized communications (check amount of overlap
        % along x and y directions)
        overlap_q = max(dims_overlap_q) - dims_q;
        u1 = comm2d_reduce(u1, overlap_q, Qyp, Qxp);     % see how to update g1 from here... 
        u1 = u1(overlap_q(1)+1:end, overlap_q(2)+1:end); % reduce to facet without overlap
        l21_cell = gop(@plus, l21_, 1);                  % reduce l21_ on the worker 1
    end

    
    %% L2 ball projection update
    % parfor can be used here, or even spmd (to be seen)
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
    g2=[]; Fx=[];
    
    %% Update primal gradient
    for i = 1 : length(v0)
        temp = v0{i};
        for j = 1 : c_chunk(i)
            g0{i}(:,:,j) = reshape(temp(:,j),M,N);
        end
    end
    
    % [P.-A.] to be completely changed here... [sub-optimal]
    g1 = cell(Qc, 1);
    for i = 1 : Qc
        g1{i} = zeros(M, N);
        for q = 1 : Q
            qi = (q-1)*Qc + i;
            g1{i}(I(q, :) + dims(q, :)) = u1{qi}; % reconstruct the full image on the master node... absolutely sub-optimal
        end
    end
    
    g = sigma00*Spt(g0) + sigma11*Spt(g1) + sigma22*Ftx; % idem, sub-optimal
    % see structure of Sp and Spt
    % Free memory
    g0=[]; g1=[]; Ftx=[];
    
    % [P.-A.]
    l21(t) = l21_cell{1};
    
    %% Display
    if ~mod(t,25)
        
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
            fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear(t), l21(t), rel_fval(t));
            fprintf(' epsilon = %e, residual = %e\n', norm(epsilon_check),norm(residual_check));
            fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
        end
        
        %diary('WB_new');
        
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
    if param.step_flag && rel_fval(t) < param.reweight_rel_obj
        reweight_steps = [t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size)];
        param.step_flag = 0;
    end
    
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
            (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
            t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);
        
        xsol_split = Sp(xsol);
        for i = 1 : length(xsol_split)
            temp = xsol_split{i};
            sol = reshape(temp(:),numel(temp(:))/c_chunk(i),c_chunk(i));
            [~,S0,~] = svd(sol,'econ');
            d_val0 = abs(diag(S0));
            weights0{i} = reweight_alpha ./ (reweight_alpha + d_val0);
            weights0{i}(d_val0 > max(d_val0) * param.reweight_abs_of_max) = 0;
        end
        
        % [P.-A.] modification for sdwt2
        % for k = 1 : P
        %     for i = 1 : length(xhat_split)
        %         d_val1 = sqrt(sum(abs(Psit{k}(xsol_split{i})).^2,2));
        %         weights1{k}{i} = reweight_alpha ./ (reweight_alpha + d_val1);
        %         weights1{k}{i}(d_val1 > max(d_val1) * param.reweight_abs_of_max) = 0;
        %     end
        % ends
        % communicate x_sol_split to the workers [sub-optimal] (ideally, x_sol lives already on the worker, update performed there...)
        for qi = 1:numworkers
            [i, q] = ind2sub([Qc, Q], qi);
            x_sol_overlap{qi} = xsol_split(I_overlap_ref(q, 1)+1:I_overlap_ref(q, 1)+dims_overlap_ref(q, 1), ...
            I_overlap_ref(q, 2)+1:I_overlap_ref(q, 2)+dims_overlap_ref(q, 2), spectral_facets{i});
        end

        spmd % encapsulate details in a function?
            zerosNum = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding (to be checked again...)
            x_ = zeros([zerosNum, size(x_overlap, 3)]);
            x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_sol_overlap;
            w = zeros(size(v1));
            for l = 1 : size(x_, 3)
                w(:, l) = sdwt2_sara(x_(:, :, l), Iq, dims_q, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q);
            end
            d_val1 = sqrt(sum(abs((w)).^2,2));
            weights1 = reweight_alphap ./ (reweight_alphap + d_val1);
            reweight_alphap = reweight_alpha_ff .* reweight_alphap;          
        end
        %-- end modifications

        reweight_alpha = reweight_alpha_ff .* reweight_alpha;
        
        if (reweight_step_count >= param.total_reweights) %|| (norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check))
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
        
        %         if (reweight_step_count == 1) || (~mod(reweight_step_count,5))
        %             save(['./reweight_res/temp_result_',num2str(reweight_step_count),'.mat'],'-v7.3', 'xsol', 'v0', 'v1', 'v2', 'g', 'weights0_old', 'weights1_old', 'weights0', 'weights1', 'proj', 't_block','reweight_alpha', 'epsilon', 't', 'rel_fval', 'nuclear', 'l21', 'norm_res', 'res');
        %         end
        
        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;
        rw_counts = rw_counts + 1;
        
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
        
    end
    
    %toc;
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

%Final log
if (param.verbose > 0)
    if (flag == 1)
        fprintf('Solution found\n');
        fprintf(' Relative variation = %e\n', rel_fval(t));
        fprintf(' Final residual = %e\n', residual_check);
        fprintf(' epsilon = %e\n', epsilon_check);
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf(' Relative variation = %e\n', rel_fval(t));
        fprintf(' Final residual = %e\n', residual_check);
        fprintf(' epsilon = %e\n', epsilon_check);
    end
end
end

% function [v1_, u1_, l21_] = run_par_waverec(v1_, Psit, Psi, xhat_split, weights1_, beta1, c)

% l21_ = 0;
% for i = 1 : length(xhat_split)
%     r1 = v1_{i} +  Psit(xhat_split{i});
%     l2 = sqrt(sum(abs(r1).^2,2));
%     l2_soft = max(l2 - beta1*weights1_{i}, 0)./ (l2+eps);
%     v1_{i} = r1 - (repmat(l2_soft,1,c(i)) .* r1);
%     u1_{i} = Psi(v1_{i});
    
%     % local L21 norm of current solution
%     l21_ = l21_ + norm(l2(:),1);
% end

end
