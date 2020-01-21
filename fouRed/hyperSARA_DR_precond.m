function [xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = hyperSARA_DR_precond(y, epsilon, A, At, H, W, pU, T, Wm, Psi, Psit, param, reduction_version, realdatablocks, fouRed_gamma)

% This function solves:
%
% min || X ||_* + lambda * ||Psit(X)||_2,1   s.t.  || Y - A(X) ||_2 <= epsilon and x>=0
%
% Author: Ming Jiang

% Useful functions for the projection
sc = @(z,eps) z*min(eps/norm(z(:)), 1); % scaling
hardt = @(z) max(real(z), 0); %thresholding negative values

% Number of nodes
R = length(H{1});
P = length(Psit);

% Number of channels
c = length(y);

% Variable flag for the case where W is present
flagW = 0;
if ~isempty(W)
    flagW = 1;
end

% number of over-sampled pixels
if flagW
    No = size(W{1}{1}, 1);
else
    No = size(H{1}{1}, 2);
end

% number of pixels
[M, N] = size(At(zeros(No, 1)));
% for i = 1 : c
%     x0(:,:,i) = reshape(X0(:,i),M,N);
% end

%Initializations.
if isfield(param,'init_xsol')
    xsol = param.init_xsol;
    fprintf('xsol uploaded \n\n')
else
    xsol = zeros(M,N,c);
    fprintf('xsol NOT uploaded \n\n')
end
sol = reshape(xsol(:),numel(xsol(:))/c,c);

%Initial dual variables
if isfield(param,'init_v0')
    v0 = param.init_v0;
    fprintf('v0 uploaded \n\n')
else
    v0 = zeros(size(sol));
    fprintf('v0 NOT uploaded \n\n')
end

if isfield(param,'init_weights0')
    weights0 = param.init_weights0;
    fprintf('weights0 uploaded \n\n')
else
    weights0 = ones(c,1);
    fprintf('weights0 NOT uploaded \n\n')
end

% %Initial dual variables
% if isfield(param,'init_v1')
%     v1 = param.init_v1;
%     fprintf('v1 uploaded \n\n')
% else
%     v1 = zeros([size(Psit(xsol(:,:,1)),1) c]);
%     fprintf('v1 NOT uploaded \n\n')
% end
% 
% if isfield(param,'weights1')
%     weights1 = param.weights1;
%     fprintf('weights1 uploaded \n\n')
% else
%     weights1 = ones(size(v1,1),1);
%     fprintf('weights1 NOT uploaded \n\n')
% end

%Initial dual variables
if isfield(param,'init_v1')
    l21_cell = cell(P, 1);
    u1 = cell(P, 1);
    v1 = param.init_v1;
    for k = 1:P
        % initial L1 descent step
        u1{k} = zeros(size(Psi{k}(v1{k})));
    end
    fprintf('v1 uploaded \n\n')
else
    l21_cell = cell(P, 1);
    u1 = cell(P, 1);
    v1 = cell(P, 1);
    for k = 1:P
        % start from zero solution
        v1{k} = zeros(size(Psit{k}(xsol)));
        
        % initial L1 descent step
        u1{k} = zeros(size(Psi{k}(v1{k})));
    end
    fprintf('v1 NOT uploaded \n\n')
end


if isfield(param,'init_weights1')
    weights1 = param.init_weights1;
    fprintf('weights1 uploaded \n\n')
else
    weights1 = cell(P, 1);
    for k = 1:P
        weights1{k} = ones(size(v1{k},1),1);
    end
    fprintf('weights1 NOT uploaded \n\n')
end


if isfield(param,'init_v2')
    v2 = param.init_v2;
    r2 = v2;
    fprintf('v2 uploaded \n\n')
else
    for i = 1 : c
        v2{i} = cell(R,1);
        for j = 1 : R
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
if isfield(param,'precondition')
    precondition = param.precondition;
else
    precondition = 0;
end

if precondition
    if isfield(param,'init_proj')
        proj = param.init_proj;
        fprintf('proj uploaded \n\n')
    else
        for i = 1 : c
            proj{i} = cell(R,1);
            Fx = A(xsol(:,:,i));
            for j = 1 : R
                if flagW
                    r2{i}{j} = T{i}{j} .* (H{i}{j} * Fx(W{i}{j}));
                else
                    r2{i}{j} = T{i}{j} .* (H{i}{j} * Fx);
                end
                [proj{i}{j}, ~] = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, r2{i}{j}, y{i}{j}, pU{i}{j}, epsilon{i}{j}, zeros(size(y{i}{j})), param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);
            end
        end
        fprintf('proj NOT uploaded \n\n')
    end
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
    for i = 1 : c
        t_block{i} = cell(R,1);
        for j = 1 : R
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

%Step size primal
tau = 0.33;

sigma00 = tau / param.nu0;
sigma11 = tau / param.nu1;
sigma22 = tau / param.nu2;


flag = 0;

beta0 = param.gamma0 * param.nu0;
beta1 = param.gamma * param.nu1;


% Gt = cell(size(T));
% for i = 1 : c
%     Gt{i} = cell(R,1);
%     for j = 1 : R
%         Gt{i}{j} = T{i}{j}';
%     end
% end

A = afclean(A);
At = afclean(At);
% Psi = afclean(Psi);
% Psit = afclean(Psit);
for k = 1 : P
    Psi{k} = afclean(Psi{k});
    Psit{k} = afclean(Psit{k});
end


% Main loop. Sequential.
%maxNumCompThreads(12);
delete(gcp('nocreate'));
util_create_pool(12);

start_loop = tic;
profile on
for t = t_start : param.max_iter
    
    start_iter = tic;
    
    %% Primal update
%     prev_xsol = xsol;
%     xsol = hardt(xsol - g);
%     xhat = 2*xsol - prev_xsol;
    [xsol, xhat, rel_x, norm_x] = update_primal(xsol, g);
    
    %% Relative change of objective function
    rel_fval(t) = sqrt(rel_x/norm_x);
    
    %% Dual variables update
    %% Nuclear norm function update
    xhatm = reshape(xhat(:),numel(xhat(:))/c,c);
    [U0,S0,V0] = svd(v0 + xhatm,'econ');
    v0 = v0 + xhatm - (U0*diag(max(diag(S0) - beta0 * weights0, 0))*V0');
%     [v0, g0] = run_par_nuclear(v0, xhat, weights0, beta0, sigma00);
    
    %% L-2,1 function update
    for k = 1:P
        f(k) = parfeval(@run_par_waverec, 3, v1{k}, Psit{k}, Psi{k}, xhat, weights1{k}, beta1,c);
    end
%     r1 = v1 +  Psit(xhat);
%     l2 = sqrt(sum(abs(r1).^2,2));
%     l2_soft = max(l2 - beta1*weights1, 0)./ (l2+eps);
%     v1 = r1 - (l2_soft .* r1);
    
    %% L2 ball projection update
    [v2, g2, proj, norm_res, norm_residual_check_c, norm_epsilon_check_c, norm_residual_check_a, norm_epsilon_check_a]...
                = update_data_fidelity_dr_block_new(v2, y, xhat, proj, A, At, H, W, T, Wm, pU, epsilon, ...
                param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps, sigma22, precondition, reduction_version, realdatablocks); % *_dr version when no blocking

    
%     counter = 1;
%     norm_residual_check_c = 0;
%     norm_epsilon_check_c = 0;
%     norm_residual_check_a = 0;
%     norm_epsilon_check_a = 0;
%     for i = 1 : c
%         Fx = A(xhat(:,:,i));
%         g2 = zeros(No,1);
%         for j = 1 : R
%             if reduction_version == 1
%                 if flagW
%                     tmp = FT2(real(At(H{i}{j} * Fx(W{i}{j}))));
%                 else
%                     tmp = FT2(real(At(H{i}{j} * Fx)));
%                 end
%                 tmp = tmp(:);
%                 r2{i}{j} = T{i}{j} .* tmp(Wm{i}{j});
%             elseif reduction_version == 2
%                 if flagW
%                     r2{i}{j} = T{i}{j} .* (H{i}{j} * Fx(W{i}{j}));
%                 else
%                     r2{i}{j} = T{i}{j} .* (H{i}{j} * Fx);
%                 end
%             end
%           
%             if precondition
%                 [proj{i}{j}, ~] = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, r2{i}{j}, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);
%                 v2{i}{j} = v2{i}{j} + pU{i}{j} .* r2{i}{j} - pU{i}{j} .* proj{i}{j};
%             else
%                 v2{i}{j} = v2{i}{j} + r2{i}{j} - y{i}{j} - sc(v2{i}{j} + r2{i}{j} - y{i}{j}, epsilon{i}{j});
%             end
%             if reduction_version == 1
%                 tmp = zeros(size(Wm{i}{j}));
%                 tmp(Wm{i}{j}) = T{i}{j} .* v2{i}{j};
%                 if flagW
%                     g2(W{i}{j}) = g2(W{i}{j}) + H{i}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
%                 else
%                     g2 = g2 + H{i}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
%                 end
%             elseif reduction_version == 2
%                 if flagW
%                     g2(W{i}{j}) = g2(W{i}{j}) + H{i}{j}' * (T{i}{j} .* v2{i}{j});
%                 else
%                     g2 = g2 + H{i}{j}' * (T{i}{j} .* v2{i}{j});
%                 end
%             end
%             
%             % norm of residual
%             norm_res{i}{j} = norm(r2{i}{j} - y{i}{j});
% %             residual_check(counter) = norm_res{i}{j};
% %             epsilon_check(counter) = epsilon{i}{j};
% %             counter = counter + 1;
%             
%             % Only for real data
%             if (realdatablocks == 2 && j == 1) || (realdatablocks == 9 && (j == 1 || j == 2))
%                 norm_residual_check_c = norm_residual_check_c + norm_res{i}{j}^2;
%                 norm_epsilon_check_c = norm_epsilon_check_c + epsilon{i}{j}^2;
%             else
%                 norm_residual_check_a = norm_residual_check_a + norm_res{i}{j}^2;
%                 norm_epsilon_check_a = norm_epsilon_check_a + epsilon{i}{j}^2;
%             end
%         end
%         Ftx(:,:,i) = real(At(g2));
%     end
%     
%     % Only for real data
%     norm_epsilon_check_c = sqrt(norm_epsilon_check_c);
%     norm_residual_check_c = sqrt(norm_residual_check_c);
%     norm_epsilon_check_a = sqrt(norm_epsilon_check_a);
%     norm_residual_check_a = sqrt(norm_residual_check_a);
%     
%     % Free memory
%     g2=[]; Fx=[];
    
    %% Update primal gradient
    g0 = reshape(v0,M,N,c); 
    g1 = zeros(size(xsol));
    for k = 1:P
        [idx, v1_, u1_, l21_] = fetchNext(f);
        v1{idx} = v1_;
        u1{idx} = u1_;
        l21_cell{idx} = l21_;
        
        g1 = g1 + u1{idx};
    end
    
    g = sigma00*g0 + sigma11*g1 + g2;
    
    end_iter(t) = toc(start_iter); 
    fprintf('Iter = %i, Time = %e\n',t,end_iter(t));
    
    %% Display
    if ~mod(t,25)
        
        %SNR
%         sol = reshape(xsol(:),numel(xsol(:))/c,c);
%         SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
%         psnrh = zeros(c,1);
%         for i = 1:c
%             psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
%         end
%         SNR_average = mean(psnrh);
        

        xhatm = reshape(xsol,numel(xsol)/c,c);
        [~,S0,~] = svd(xhatm,'econ');
        nuclear = norm(diag(S0),1);
        
        l21 = sum(cell2mat(l21_cell));
%         l2 = sqrt(sum(abs(Psit(xsol)).^2,2));
%         l21 = norm(l2(:),1);
        
        %Log
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
            fprintf('epsilon_c = %e, residual_c = %e\n', norm_epsilon_check_c, norm_residual_check_c);
            fprintf('epsilon_a = %e, residual_a = %e\n', norm_epsilon_check_a, norm_residual_check_a);
        end
        
%           diary('WB_new');
        fitswrite(xsol, ['hyper_xsol_it', num2str(t), '_gamma', num2str(param.gamma), '_', num2str(realdatablocks),...
            'b_fouRed', num2str(reduction_version), '_th', num2str(fouRed_gamma), '.fits']);
    end
    
    %% Global stopping criteria
    %     prod(prod(residual_check < param.adapt_eps_tol_out*epsilon_check)) && prod(prod(residual_check > param.adapt_eps_tol_in*epsilon_check))
    if t>1 && rel_fval(t) < param.rel_obj && reweight_step_count > param.total_reweights && ...
            (norm_residual_check_c <= param.adapt_eps_tol_out*norm_epsilon_check_c) && ...
            (norm_residual_check_a <= param.adapt_eps_tol_out*norm_epsilon_check_a)
        flag = 1;
        break;
    end
    
    %% Update epsilons
    if param.use_adapt_eps && t > param.adapt_eps_start
        for i = 1 : c
            for  j = 1 : R
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
%     if param.step_flag && rel_fval(t) < param.reweight_rel_obj % && (norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check))
%         reweight_steps = t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size);
%         param.step_flag = 0;
%     end
    
    % Only for real data
    if (param.step_flag && rel_fval(t) < param.reweight_rel_obj && ...
            norm_residual_check_c <= param.adapt_eps_tol_out*norm_epsilon_check_c && norm_residual_check_a <= param.adapt_eps_tol_out*norm_epsilon_check_a && t > 300)
        reweight_steps = (t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size));
        param.step_flag = 0;
    end
    
%     if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
%             (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
%             t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
    
    % Only for real data
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
          (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
          t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr && ...
          norm_residual_check_c <= param.adapt_eps_tol_out*norm_epsilon_check_c && ...
          norm_residual_check_a <= param.adapt_eps_tol_out*norm_epsilon_check_a) || ...
          (param.use_reweight_eps && t - reweight_last_step_iter > param.rw_tol)
        %(norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check)) && ...
        
%         weights0_old = weights0;
%         weights1_old = weights1;
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);
        
        sol = reshape(xsol(:),numel(xsol(:))/c,c);
        [~,S0,~] = svd(sol,'econ');
        d_val0 = abs(diag(S0));
        weights0 = reweight_alpha ./ (reweight_alpha + d_val0);
        weights0(d_val0 > max(d_val0) * param.reweight_abs_of_max) = 0;
        
%         d_val1 = sqrt(sum(abs(Psit(xsol)).^2,2));
%         weights1 = reweight_alpha ./ (reweight_alpha + d_val1);
%         weights1(d_val1 > max(d_val1) * param.reweight_abs_of_max) = 0;
        
        for k = 1 : P
            d_val1 = sqrt(sum(abs(Psit{k}(xsol)).^2,2));
            weights1{k} = reweight_alpha ./ (reweight_alpha + d_val1);
            weights1{k}(d_val1 > max(d_val1) * param.reweight_abs_of_max) = 0;
        end
        reweight_alpha = reweight_alpha_ff .* reweight_alpha;
        
        if (reweight_step_count >= param.total_reweights) %|| (norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check))
            param.reweight_max_reweight_itr = t+1;
            fprintf('\n\n No more reweights \n\n');
            %break;
        end
        
        % Calculate residual images:
        for i = 1 : c
            Fx = A(xsol(:,:,i));
            g2 = zeros(No,1);
            for j = 1 : R
                if reduction_version == 1
                    if flagW
                        tmp = FT2(real(At(H{i}{j} * Fx(W{i}{j}))));
                    else
                        tmp = FT2(real(At(H{i}{j} * Fx)));
                    end
                    tmp = tmp(:);
                    res_f{i}{j} = y{i}{j} - T{i}{j} .* tmp(Wm{i}{j});
                    tmp = zeros(size(Wm{i}{j}));
                    tmp(Wm{i}{j}) = T{i}{j} .* res_f{i}{j};
                    if flagW
                        g2(W{i}{j}) = g2(W{i}{j}) + H{i}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
                    else
                        g2 = g2 + H{i}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
                    end
                elseif reduction_version == 2
                    if flagW
                        res_f{i}{j} = y{i}{j} - T{i}{j} .* (H{i}{j} * Fx(W{i}{j}));
                        g2(W{i}{j}) = g2(W{i}{j}) + H{i}{j}' * (T{i}{j} .* res_f{i}{j});
                    else
                        res_f{i}{j} = y{i}{j} - T{i}{j} .* (H{i}{j} * Fx);
                        g2 = g2 + H{i}{j}' * (T{i}{j} .* res_f{i}{j});
                    end
                end
            end
%             res(:,:,i) = real(At(g2));
        end
        
%         fprintf('\n\n\n\n\n\n\n Performed reweight no %d \n\n\n\n\n', reweight_step_count);
        fitswrite(xsol, ['hyper_xsol_it', num2str(t), '_reweight', num2str(reweight_step_count), '_gamma', num2str(param.gamma)...
            '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_th', num2str(fouRed_gamma), '.fits']);
        
        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;
        rw_counts = rw_counts + 1;
        
    end
    
    %toc;
end
profile off
end_loop = toc(start_loop)
profsave(profile('info'),'HyperSARA_DR_profile_results')

% Calculate residual images:
for i = 1 : c
    Fx = A(xsol(:,:,i));
    g2 = zeros(No,1);
    for j = 1 : R
        if reduction_version == 1
            if flagW
                tmp = FT2(real(At(H{i}{j} * Fx(W{i}{j}))));
            else
                tmp = FT2(real(At(H{i}{j} * Fx)));
            end
            tmp = tmp(:);
            res_f{i}{j} = y{i}{j} - T{i}{j} .* tmp(Wm{i}{j});
            tmp = zeros(size(Wm{i}{j}));
            tmp(Wm{i}{j}) = T{i}{j} .* res_f{i}{j};
            if flagW
                g2(W{i}{j}) = g2(W{i}{j}) + H{i}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
            else
                g2 = g2 + H{i}{j} * A(real(IFT2(reshape(tmp, Ny, Nx))));
            end
        elseif reduction_version == 2
            if flagW
                res_f{i}{j} = y{i}{j} - T{i}{j} .* (H{i}{j} * Fx(W{i}{j}));
                g2(W{i}{j}) = g2(W{i}{j}) + H{i}{j}' * (T{i}{j} .* res_f{i}{j});
            else
                res_f{i}{j} = y{i}{j} - T{i}{j} .* (H{i}{j} * Fx);
                g2 = g2 + H{i}{j}' * (T{i}{j} .* res_f{i}{j});
            end
        end
    end
    res(:,:,i) = real(At(g2));
end

%Final log

xhatm = reshape(xsol,numel(xsol)/c,c);
[~,S0,~] = svd(xhatm,'econ');
nuclear = norm(diag(S0),1);

l21 = sum(cell2mat(l21_cell));

% %SNR
% sol = reshape(xsol(:),numel(xsol(:))/c,c);
% SNR = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)));
% psnrh = zeros(c,1);
% for i = 1:c
%     psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
% end
% SNR_average = mean(psnrh);

if (param.verbose > 0)
    if (flag == 1)
        fprintf('Solution found\n');
        fprintf('Iter %i\n',t);
        fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
        fprintf('epsilon_c = %e, residual_c = %e\n', norm_epsilon_check_c, norm_residual_check_c);
        fprintf('epsilon_a = %e, residual_a = %e\n', norm_epsilon_check_a, norm_residual_check_a);
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf('Iter %i\n',t);
        fprintf('N-norm = %e, L21-norm = %e, rel_fval = %e\n', nuclear, l21, rel_fval(t));
        fprintf('epsilon_c = %e, residual_c = %e\n', norm_epsilon_check_c, norm_residual_check_c);
        fprintf('epsilon_a = %e, residual_a = %e\n', norm_epsilon_check_a, norm_residual_check_a);
    end
end

end

function [v1_, u1_, l21_] = run_par_waverec(v1_, Psit, Psi, xhat, weights1_, beta1, c)

r1 = v1_ +  Psit(xhat);
l2 = sqrt(sum(abs(r1).^2,2));
l2_soft = max(l2 - beta1*weights1_, 0)./ (l2+eps);
v1_ = r1 - (repmat(l2_soft,1,c) .* r1);
u1_ = Psi(v1_);

% local L21 norm of current solution
l21_ = norm(l2(:),1);
end
