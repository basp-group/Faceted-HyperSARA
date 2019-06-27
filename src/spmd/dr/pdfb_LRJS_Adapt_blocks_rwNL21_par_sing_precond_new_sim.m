function [xsol,v0,v1,v2,g,weights0,weights1,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res] = pdfb_LRJS_Adapt_blocks_rwNL21_par_sing_precond_new_sim(y, imsize, epsilon, A, At, pU, T, W, Psi, Psit, param,X0)

% This function solves:
%
% min || X ||_* + lambda * ||Psit(X)||_2,1   s.t.  || Y - A(X) ||_2 <= epsilon and x>=0
%
% Author: Ming Jiang

% Useful functions for the projection
sc = @(z,eps) z*min(eps/norm(z(:)), 1); % scaling
hardt = @(z) max(real(z), 0); %thresholding negative values

% Number of nodes
R = length(T{1});
P = length(Psit);

% number of pixels
M = imsize(1);
N = imsize(2);
c = imsize(3);

for i = 1 : c
    x0(:,:,i) = reshape(X0(:,i),M,N);
end

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

%Initial dual variables
if isfield(param,'init_v1')
    l2 = cell(P,1);
    l21_cell = cell(P, 1);
    u1 = cell(P, 1);
    v1 = cell(P, 1);
    v1 = param.init_v1;
    for k = 1:P
        % initial L1 descent step
        u1{k} = zeros(size(Psi{k}(v1{k})));
    end
    fprintf('v1 uploaded \n\n')
else
    l2 = cell(P,1);
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
    for i = 1 : c
        u2{i} = cell(R, 1);
        for j = 1 : R
            u2{i}{j} = zeros(size(T{i}{j}, 2), 1);
        end
    end
    r2 = v2;
    fprintf('v2 uploaded \n\n')
else
    for i = 1 : c
        v2{i} = cell(R,1);
        u2{i} = cell(R,1);
        for j = 1 : R
            v2{i}{j} = zeros(length(y{i}{j}) ,1);
            u2{i}{j} = zeros(size(T{i}{j}, 2), 1);
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

g0 = zeros(size(xsol));
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
            Fx = A{i}(xsol(:,:,i));
            for j = 1 : R
                r2{i}{j} = T{i}{j} .* Fx(W{i}{j});
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


% Gt = cell(size(T));
% for i = 1 : c
%     Gt{i} = cell(R,1);
%     for j = 1 : R
%         Gt{i}{j} = T{i}{j}';
%     end
% end

% A = afclean(A);
% At = afclean(At);
for k = 1 : P
    Psi{k} = afclean(Psi{k});
    Psit{k} = afclean(Psit{k});
end


% Main loop. Sequential.
%maxNumCompThreads(12);
util_create_pool(24);

for t = t_start : param.max_iter
    
    %fprintf('Iter %i\n',t);
    %tic;
    
    %% Primal update
    prev_xsol = xsol;
    xsol = hardt(xsol - g);
    xhat = 2*xsol - prev_xsol;
    
    %% Relative change of objective function
    rel_fval(t) = norm(xsol(:) - prev_xsol(:))/norm(xsol(:));
    % Free memory
    prev_xsol = [];
    
    %% Dual variables update
    %% Nuclear norm function update
    xhatm = reshape(xhat(:),numel(xhat(:))/c,c);
    [U0,S0,V0] = svd(v0 + xhatm,'econ');
    nuclear(t) = norm(diag(S0),1);
    v0 = v0 + xhatm - (U0*diag(max(diag(S0) - beta0 * weights0, 0))*V0');
    % Free memory
    U0=[]; S0=[]; V0=[]; xhatm = [];
    
    
    %% L-2,1 function update
    for k = 1:P
        f(k) = parfeval(@run_par_waverec, 3, v1{k}, Psit{k}, Psi{k}, xhat, weights1{k}, beta1,c);
    end
    
    %% L2 ball projection update
    counter = 1;
    for i = 1 : c
        Fx = A{i}(xhat(:,:,i));
        g2 = zeros(M * N,1);
        for j = 1 : R
            r2{i}{j} = T{i}{j} .* Fx(W{i}{j});
            if precondition
                [proj{i}{j}, ~] = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, r2{i}{j}, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);
                v2{i}{j} = v2{i}{j} + pU{i}{j} .* r2{i}{j} - pU{i}{j} .* proj{i}{j};
            else
                v2{i}{j} = v2{i}{j} + r2{i}{j} - y{i}{j} - sc(v2{i}{j} + r2{i}{j} - y{i}{j}, epsilon{i}{j});
            end
                
            u2{i}{j} = T{i}{j} .* v2{i}{j};
            g2(W{i}{j}) = g2(W{i}{j}) + u2{i}{j};
            
            % norm of residual
            norm_res{i}{j} = norm(r2{i}{j} - y{i}{j});
            residual_check(counter) = norm_res{i}{j};
            epsilon_check(counter) = epsilon{i}{j};
            counter = counter + 1;
        end
        Ftx(:,:,i) = real(At{i}(g2));
    end
    
    % Free memory
    g2=[]; Fx=[];
    
    %% Update primal gradient
    for i = 1 : c
        g0(:,:,i) = reshape(v0(:,i),M,N);
    end
    
    g1 = zeros(size(xsol));
    for k = 1:P
        [idx, v1_, u1_, l21_] = fetchNext(f);
        v1{idx} = v1_;
        u1{idx} = u1_;
        l21_cell{idx} = l21_;
        
        g1 = g1 + u1{idx};
    end
    
    g = sigma00*g0 + sigma11*g1 + sigma22*Ftx;
    % Free memory
    g0=[]; g1=[]; Ftx=[];
    
    l21(t) = sum(cell2mat(l21_cell));
    
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
        
          diary('WB_new');
    end
    
    %% Save
    if ~mod(t,200000)
        
        % Calculate residual images:
        for i = 1 : c
            Fx = A{i}(xsol(:,:,i));
            g2 = zeros(N * N,1);
            for j = 1 : R
                res_f{i}{j} = y{i}{j} - T{i}{j} .* Fx(W{i}{j});
                u2{i}{j} = T{i}{j} .* res_f{i}{j};
                g2(W{i}{j}) = g2(W{i}{j}) + u2{i}{j};
            end
            res(:,:,i) = real(At{i}(g2));
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
    %if (param.use_reweight_steps && t == param.reweight_steps(reweight_step_count)) || ...
    %        (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
    %        prod(prod(residual_check < param.adapt_eps_tol_out*epsilon_check)) && prod(prod(residual_check > param.adapt_eps_tol_in*epsilon_check))  && ...
    %       t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
    
    if param.step_flag && rel_fval(t) < param.reweight_rel_obj % && (norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check))
        reweight_steps = [t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size)];
        param.step_flag = 0;
    end
    
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
            (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
            t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
        %(norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check)) && ...
        
        weights0_old = weights0;
        weights1_old = weights1;
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);
        
        sol = reshape(xsol(:),numel(xsol(:))/c,c);
        [~,S0,~] = svd(sol,'econ');
        d_val0 = abs(diag(S0));
        weights0 = reweight_alpha ./ (reweight_alpha + d_val0);
        weights0(d_val0 > max(d_val0) * param.reweight_abs_of_max) = 0;
        
        for k = 1 : P
            d_val1 = sqrt(sum(abs(Psit{k}(xsol)).^2,2));
            weights1{k} = reweight_alpha ./ (reweight_alpha + d_val1);
            weights1{k}(d_val1 > max(d_val1) * param.reweight_abs_of_max) = 0;
        end
        reweight_alpha = reweight_alpha_ff .* reweight_alpha;
        
        if (reweight_step_count >= param.total_reweights) %|| (norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check))
            param.reweight_max_reweight_itr = t+1;
            fprintf('\n\n No more reweights \n\n');
            break;
        end
        
        % Calculate residual images:
        for i = 1 : c
            Fx = A{i}(xsol(:,:,i));
            g2 = zeros(M * N,1);
            for j = 1 : R
                res_f{i}{j} = y{i}{j} - T{i}{j} .* Fx(W{i}{j});
                u2{i}{j} = T{i}{j} .* res_f{i}{j};
                g2(W{i}{j}) = g2(W{i}{j}) + u2{i}{j};
            end
            res(:,:,i) = real(At{i}(g2));
        end
        
        %         if (reweight_step_count == 1) || (~mod(reweight_step_count,5))
        %             save(['./reweight_res/temp_result_',num2str(reweight_step_count),'.mat'],'-v7.3', 'xsol', 'v0', 'v1', 'v2', 'g', 'weights0_old', 'weights1_old', 'weights0', 'weights1', 'proj', 't_block','reweight_alpha', 'epsilon', 't', 'rel_fval', 'nuclear', 'l21', 'norm_res', 'res');
        %         end
        
        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;
        rw_counts = rw_counts + 1;
        
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
        
    end
    
    %toc;
end

% Calculate residual images:
for i = 1 : c
    Fx = A{i}(xsol(:,:,i));
    g2 = zeros(M * N,1);
    for j = 1 : R
        res_f{i}{j} = y{i}{j} - T{i}{j} .* Fx(W{i}{j});
        u2{i}{j} = T{i}{j} .* res_f{i}{j};
        g2(W{i}{j}) = g2(W{i}{j}) + u2{i}{j};
    end
    res(:,:,i) = real(At{i}(g2));
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

function [v1_, u1_, l21_] = run_par_waverec(v1_, Psit, Psi, xhat, weights1_, beta1, c)

r1 = v1_ +  Psit(xhat);
l2 = sqrt(sum(abs(r1).^2,2));
l2_soft = max(l2 - beta1*weights1_, 0)./ (l2+eps);
v1_ = r1 - (repmat(l2_soft,1,c) .* r1);
u1_ = Psi(v1_);

% local L21 norm of current solution
l21_ = norm(l2(:),1);
end