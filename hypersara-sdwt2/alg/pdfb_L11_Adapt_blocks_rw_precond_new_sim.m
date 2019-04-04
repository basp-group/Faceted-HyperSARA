function [xsol,v1,v2,g,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,l11,norm_res,res,end_iter] = pdfb_L11_Adapt_blocks_rw_precond_new_sim(y, epsilon, A, At, pU, G, W, Psi, Psit, param, X0, nu2)

% This function solves:
%
%   min lambda * ||Psit(X)||_2,1   s.t.  || Y - A(X) ||_2 <= epsilon and x>=0
%
% Author: Abdullah Abdulaziz

% Useful functions for the projection
% sc = @(z,eps) z*min(eps/norm(z(:)), 1); % scaling
hardt = @(z) max(real(z), 0); %thresholding negative values

c = size(y,2);
P = length(Psit);

% oversampling vectorized data length
No = size(W{1}{1}, 1);

% number of pixels
[M, N] = size(At(zeros(No, 1)));

x0 = reshape(X0,M,N);

param.nu2 = nu2;

%Initializations.
if isfield(param,'init_xsol')
    xsol = param.init_xsol;
    fprintf('xsol uploaded \n\n')
else
    xsol = zeros(M,N,c);
    fprintf('xsol NOT uploaded \n\n')
end

%Initial dual variables
if isfield(param,'init_v1')
    v1 = param.init_v1;
    fprintf('v1 uploaded \n\n')
else
    v1 = zeros([size(Psit(xsol(:,:,1)),1) c]);
    fprintf('v1 NOT uploaded \n\n')
end

if isfield(param,'weights1')
    weights1 = param.weights1;
    fprintf('weights1 uploaded \n\n')
else
    weights1 = ones(size(v1,1),1);
    fprintf('weights1 NOT uploaded \n\n')
end


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
    rw_counts = 1;
    fprintf('t t_block uploaded \n\n')
else
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
sigma1 = 1.0/param.nu1;
sigma2 = 1.0/param.nu2;

%Step size primal
tau = 0.99/(sigma1*param.nu1 + sigma2*param.nu2);

sigma11 = tau*sigma1;
sigma22 = tau*sigma2;

flag = 0;

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
Psi = afclean(Psi);
Psit = afclean(Psit);




% Main loop. Sequential.

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
    
    %% Dual variables update
    
    %% L-1,1 function update
    r1 = v1 + Psit(xhat);
    v1 = r1 - sign(r1).*max(abs(r1) - beta1*weights1,0);
    
    % local L11 norm of current solution
    l11(t) = norm(r1(:),1);
    
    %% L2 ball projection update
    counter = 1;
    for i = 1 : c
        Fx = A(xhat(:,:,i));
        g2 = zeros(No,1);
        for j = 1 : length(G{i})
            r2{i}{j} = G{i}{j} * Fx(W{i}{j});
            
            [proj{i}{j}, ~] = solver_proj_elipse_fb(1 ./ pU{i}{j} .* v2{i}{j}, r2{i}{j}, y{i}{j}, pU{i}{j}, epsilon{i}{j}, proj{i}{j}, param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);
            v2{i}{j} = v2{i}{j} + pU{i}{j} .* r2{i}{j} - pU{i}{j} .* proj{i}{j};
            
            % projection onto the l2-ball
            %v2{i}{j} = v2{i}{j} + r2{i}{j} - proj_l2ball(v2{i}{j} + r2{i}{j}, epsilon{i}{j}, y{i}{j});
            
            
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
    
    g = sigma11*Psi(v1) + sigma22*Ftx;
    % Free memory
    %Ftx=[]; g1=[];
    
    end_iter(t) = toc(start_iter);
    fprintf('Iter = %i, Time = %e\n',t,end_iter(t));
    
    %% Display
    if ~mod(t,25000)
        
        %SNR
        SNR = 20*log10(norm(X0(:))/norm(X0(:)-xsol(:)));
        
        %Log
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf('\n l11-norm = %e, rel_fval = %e\n\n', l11(t), rel_fval(t));
            fprintf(' epsilon = %e, residual = %e\n\n', norm(epsilon_check),norm(residual_check));
            fprintf(' SNR = %e\n\n', SNR);
            
            figure(1),
            
            subplot(1,2,1);
            imagesc(log10(max(flip(x0),0))); hold on;
            colorbar;
            axis image;
            axis off;
            colormap(cubehelix);
            caxis([-3.5, 0]);
            xlabel(['x0 - it = ',num2str(t)])
            
            subplot(1,2,2);
            imagesc(log10(max(flip(xsol),0))); hold on;
            colorbar;
            axis image;
            axis off;
            colormap(cubehelix);
            caxis([-3.5, 0]);
            xlabel(['xsol - it = ',num2str(t)])
            
            pause(0.1)
        end
        
        %diary('erw_a06_L1_5e5_6b');
    end
    
    %% Global stopping criteria
    %     if (t>1 && rel_fval(t) < param.rel_obj && reweight_step_count >= param.total_reweights && ...
    %             prod(prod(residual_check < param.adapt_eps_tol_out*epsilon_check)) && prod(prod(residual_check > param.adapt_eps_tol_in*epsilon_check)))
    %         flag = 1;
    %         break;
    %     end
    
    if (t>1 && rel_fval(t) < param.rel_obj && reweight_step_count >= param.total_reweights && ...
            (norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check)))
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
    %if (param.use_reweight_steps && t == param.reweight_steps(reweight_step_count)) || ...
    %        (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
    %        prod(prod(residual_check < param.adapt_eps_tol_out*epsilon_check)) && prod(prod(residual_check > param.adapt_eps_tol_in*epsilon_check))  && ...
    %       t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
    
    if (param.step_flag && rel_fval(t) < param.reweight_rel_obj) % && (norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check))
        reweight_steps = [t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size)];
        param.step_flag = 0;
    end
    
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
            (param.use_reweight_eps && rel_fval(t) < param.reweight_rel_obj && ...
            t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
        %(norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check)) && ...
        
        weights1_old = weights1;
        
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);
 
        
        d_val1 = abs(Psit(xsol));
        weights1 = reweight_alpha ./ (reweight_alpha + d_val1);
        weights1(d_val1 > max(d_val1) * param.reweight_abs_of_max) = 0;
        reweight_alpha = reweight_alpha_ff .* reweight_alpha;
        
        
        
        if reweight_step_count > param.total_reweights
            param.reweight_max_reweight_itr = t+1;
            fprintf('\n\n No more reweights \n\n');
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
        
        %fitswrite(xsol,['./reweight_res/rw_sol_',num2str(reweight_step_count),'.fits']);
        %fitswrite(res,['./reweight_res/rw_res_',num2str(reweight_step_count),'.fits']);
        %         if (reweight_step_count == 30) || (reweight_step_count == 35) || (reweight_step_count == 40) || ...
        %            (reweight_step_count == 45) || (reweight_step_count == 50) %(~mod(reweight_step_count,5))
        %             save(['./reweight_res_l1/temp_result_' num2str(ch) '_' num2str(reweight_step_count) '.mat'],'-v7.3', 'xsol', 'v1', 'v2', 'g', 'weights1_old', 'weights1', 'proj', 't_block', 'reweight_alpha', 'epsilon', 'rel_fval', 'l11', 'norm_res', 'res');
        %         end
        
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
