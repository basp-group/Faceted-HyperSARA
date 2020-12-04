function [xsol,param,v1,v2,g,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_val,l11,norm_res,res,t_l11,t_master,end_iter] = ...
    sara(y, epsilon, A, At, pU, G, W, Psi, Psit, param, ch, init_file_name, name, x0)

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


%Initializations.
init_flag = isfile(init_file_name)
if init_flag
    init_m = matfile(init_file_name)
end

if init_flag
    xsol = init_m.xsol;
    param = init_m.param;
    epsilon = init_m.epsilon;
    g = init_m.g;
    fprintf('xsol, g, param and epsilon uploaded \n\n')
else
    xsol = zeros(M,N,c);
    g = zeros(size(xsol));
    %! [P.-A.]
    % xsol = param.xsol_Arwa; 
    % param.xsol_Arwa = [];
    % fprintf('xsol initialized \n\n')
    fprintf('xsol and g initialized \n\n')
end

% Initial dual variables
if init_flag
    l11_cell = cell(P, 1);
    v1 = init_m.v1;
    u1 = cell(P, 1);
    for k = 1:P
        % initial L1 descent step
        u1{k} = zeros(size(Psi{k}(v1{k})));
    end
    weights1 = init_m.weights1;
    fprintf('v1, weights1 uploaded \n\n')
else
    l11_cell = cell(P, 1);
    u1 = cell(P, 1);
    v1 = cell(P, 1);
    weights1 = cell(P, 1);
    for k = 1:P
        % start from zero solution
        v1{k} = zeros(size(Psit{k}(xsol)));
        % initial L1 descent step
        u1{k} = zeros(size(Psi{k}(v1{k})));
        weights1{k} = ones(size(v1{k},1),1);
    end
    fprintf('v1, weights1 initialized \n\n')
end

if init_flag
    v2 = init_m.v2;
    for i = 1 : c
        u2{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            u2{i}{j} = zeros(size(G{i}{j}, 2), 1);
        end
    end
    r2 = v2;
    norm_res = init_m.norm_res;
    fprintf('v2, norm_res uploaded \n\n')
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
    fprintf('v2 initialized \n\n')
end


% Initialise projection
if init_flag
    proj = init_m.proj;
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
    fprintf('proj initialized \n\n')
end


if init_flag
    t_block = init_m.t_block;
    t_start = param.init_t_start;
    reweight_last_step_iter = param.init_reweight_last_iter_step;
    reweight_step_count = param.init_reweight_step_count;
    rw_counts = 1;
    fprintf('t_start, t_block uploaded \n\n')
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
    fprintf('t_start, t_block initialized \n\n')
end

count_eps_update_down = 0;
count_eps_update_up = 0;

reweight_alpha = param.reweight_alpha;
reweight_alpha_ff = param.reweight_alpha_ff;
reweight_steps = param.reweight_steps;

%%%%% ABD
%! [P.-A.] not sure this is useful
for k = 1 : P
    d_val = abs(Psit{k}(xsol));
    weights1{k} = reweight_alpha ./ (reweight_alpha + d_val);
    % weights1{k}(d_val > max(d_val) * param.reweight_abs_of_max) = 0;
end

if init_flag
    rel_val = init_m.rel_val;
    l11 = init_m.l11;
    end_iter = init_m.end_iter;
    t_master = init_m.t_master;
    t_l11 = init_m.t_l11;
    t_data = init_m.t_data;
    fprintf('rel_val, l11, end_iter, t_master, t_l11, and t_data uploaded \n\n')
else
    rel_val = zeros(param.max_iter, 1);
    l11 = zeros(param.max_iter, 1);
    end_iter = zeros(param.max_iter, 1);
    t_master = zeros(param.max_iter, 1);
    t_l11 = zeros(param.max_iter, 1);
    t_data = zeros(param.max_iter, 1);
    fprintf('rel_val, l11, end_iter, t_master, t_l11, and t_data initialized \n\n')
end

if init_flag
    SNR = 20*log10(norm(x0(:))/norm(x0(:)-xsol(:)));  
    norm_residual_check = 0;
    norm_epsilon_check = 0;
    for i = 1 : c
        for j = 1 : length(G{i})
            norm_residual_check = norm_residual_check + norm_res{i}{j}^2;
            norm_epsilon_check = norm_epsilon_check + power(epsilon{i}{j},2);
        end
    end
    norm_residual_check = sqrt(norm_residual_check);
    norm_epsilon_check = sqrt(norm_epsilon_check);
    % Log
    if (param.verbose >= 1)
        fprintf('Iter %i\n',t_start-1);
        fprintf('l11-norm = %e, rel_val = %e\n', l11(t_start-1), rel_val(t_start-1));
        fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
        fprintf(' SNR = %e\n', SNR);
    end
end

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
for k = 1 : P
    Psi{k} = afclean(Psi{k});
    Psit{k} = afclean(Psit{k});
end
Ftx = zeros(size(xsol));

% Main loop. Sequential.
%maxNumCompThreads(12);
% util_create_pool(15); %! ask Abdullah here
spmd
    dwtmode('zpd')
end

for t = t_start : param.max_iter
    
    % fprintf('Iter %i\n',t);
    start_iter = tic;
    
    %% Primal update
    tw = tic;
    prev_xsol = xsol;
    xsol = hardt(xsol - g);
    xhat = 2*xsol - prev_xsol;
    t_master(t) = toc(tw);
    
    %% Relative change of objective function
    rel_val(t) = norm(xsol(:) - prev_xsol(:))/norm(xsol(:));
    % Free memory
    prev_xsol = [];
    
    %% Dual variables update
    
    %% L-1,1 function update
    for k = 1:P
        f(k) = parfeval(@run_par_waverec, 4, v1{k}, Psit{k}, Psi{k}, xhat, weights1{k}, beta1);
    end
    
    %% L2 ball projection update
    %! [P.-A.]
    % residual_check_c = 0;
    % epsilon_check_c = 0;
    % residual_check_a = 0;
    % epsilon_check_a = 0;
    norm_residual_check = 0;
    norm_epsilon_check = 0;
    counter = 1;
    tw = tic;
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

            %! [P.-A.] what's the difference between residual_check_c and residual_check_a ??
            % if j == 1
            %     residual_check_c = residual_check_c + norm_res{i}{j}^2;
            %     epsilon_check_c = epsilon_check_c + power(epsilon{i}{j},2);
            % else
            %     residual_check_a = residual_check_a + norm_res{i}{j}^2;
            %     epsilon_check_a = epsilon_check_a + power(epsilon{i}{j},2);
            % end
            norm_residual_check = norm_residual_check + norm_res{i}{j}^2;
            norm_epsilon_check = norm_epsilon_check + power(epsilon{i}{j},2);

        end
        Ftx(:,:,i) = real(At(g2));
    end
    t_data(t) = toc(tw);
    % Free memory
    g2=[]; Fx=[];
   
    % epsilon_check_c = sqrt(epsilon_check_c);
    % residual_check_c = sqrt(residual_check_c);
    % epsilon_check_a = sqrt(epsilon_check_a);
    % residual_check_a = sqrt(residual_check_a);
    norm_epsilon_check = sqrt(norm_epsilon_check);
    norm_residual_check = sqrt(norm_residual_check);
 
    %% Update primal gradient
    g1 = zeros(size(xsol));
    for k = 1:P
        [idx, v1_, u1_, l11_, t_l11_] = fetchNext(f);
        v1{idx} = v1_;
        u1{idx} = u1_;
        l11_cell{idx} = l11_;
        t_l11(t) = t_l11(t) + t_l11_;
        g1 = g1 + u1{idx};
    end
    g = sigma11*g1 + sigma22*Ftx;
    end_iter(t) = toc(start_iter);
    t_l11(t) = t_l11(t)/P; % average compute time for the dual variable
    l11(t) = sum(cell2mat(l11_cell));
    fprintf('Iter = %i, Time = %e, t_l11 = %e, t_data = %e\n',t,end_iter(t),t_data(t),t_l11(t));

    % Free memory
    Ftx=[]; g1=[];
    
    %% Display
    if ~mod(t,100)

        % SNR
        SNR = 20*log10(norm(x0(:))/norm(x0(:)-xsol(:)));
        
        % Log
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf('l11-norm = %e, rel_val = %e\n', l11(t), rel_val(t));
            %! [P.-A.] not sure why we have residual_c / residual_a, ...
            % fprintf(' epsilon = %e, residual = %e\n', norm(epsilon_check), norm(residual_check));
            % fprintf(' epsilon_c = %e, residual_c = %e\n', epsilon_check_c, residual_check_c);
            % fprintf(' epsilon_a = %e, residual_a = %e\n', epsilon_check_a, residual_check_a);
            fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
            fprintf(' SNR = %e\n', SNR);

            for i = 1 : length(epsilon_check)
               fprintf(['eps_b' num2str(i) '= %e, res_b' num2str(i) '= %e\n'], epsilon_check(i), residual_check(i));
            end
        end
    end
    
    %% Global stopping criteria
    %     if (t>1 && rel_val(t) < param.rel_obj && reweight_step_count >= param.total_reweights && ...
    %             prod(prod(residual_check < param.adapt_eps_tol_out*epsilon_check)) && prod(prod(residual_check > param.adapt_eps_tol_in*epsilon_check)))
    %         flag = 1;
    %         break;
    %     end
    
    %! copy-paste condition from hypersara
    if t>1 && rel_val(t) < param.rel_var && reweight_step_count > param.total_reweights && ...
        (norm_residual_check <= param.adapt_eps_tol_out*norm_epsilon_check)
        flag = 1;
        break;
    end

    % %! [P.-A.] version active before I touched the code, different from hypersara. Is it normal?
    % if (t>1 && rel_val(t) < param.rel_obj && reweight_step_count > param.total_reweights && ...
    %         residual_check_c <= param.adapt_eps_tol_out*epsilon_check_c && residual_check_a <= param.adapt_eps_tol_out*epsilon_check_a)
    %     flag = 1;
    %     break;
    % end
    
    %% Update epsilons
    %! [P.-A.]
    % if param.use_adapt_eps && t > param.adapt_eps_start
    if param.use_adapt_eps && (t > param.adapt_eps_start) && (rel_val(t) < param.adapt_eps_rel_var)
        for i = 1 : c
            for  j = 1 : length(G{i})
                if  norm_res{i}{j} < param.adapt_eps_tol_in * epsilon{i}{j}
                    if t > t_block{i}{j} + param.adapt_eps_steps && rel_val(t) < param.adapt_eps_rel_var
                        epsilon{i}{j} = norm_res{i}{j} + (-norm_res{i}{j} + epsilon{i}{j}) * (1 - param.adapt_eps_change_percentage);
                        t_block{i}{j} = t;
                        count_eps_update_down = count_eps_update_down + 1;
                        fprintf('Updated  epsilon DOWN: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j},norm_res{i}{j},j,i);
                    end
                end
                
                if  norm_res{i}{j} > param.adapt_eps_tol_out * epsilon{i}{j}
                    if t > t_block{i}{j} + param.adapt_eps_steps && rel_val(t) < param.adapt_eps_rel_var
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
    %        (param.use_reweight_eps && rel_val(t) < param.reweight_rel_var && ...
    %        prod(prod(residual_check < param.adapt_eps_tol_out*epsilon_check)) && prod(prod(residual_check > param.adapt_eps_tol_in*epsilon_check))  && ...
    %       t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr)
    
    %! [P.-A.] is reweight_steps used at all? replacing by version used in hypersara
    % if (param.step_flag && (norm(residual_check) <= param.adapt_eps_tol_out*norm(epsilon_check)) && ...
    %    rel_val(t) < param.reweight_rel_var && t > 300)
    %     reweight_steps = [t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size)];
    %     param.step_flag = 0;
    % end
    if (param.step_flag && t>500) % rel_fval(t) < param.reweight_rel_var)
        reweight_steps = (t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size));
        param.step_flag = 0;
    end
    
    %! [P.-A.]
    % if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
    %         (param.use_reweight_eps && rel_val(t) < param.reweight_rel_var && ...
    %         residual_check_c <= param.adapt_eps_tol_out*epsilon_check_c && residual_check_a <= param.adapt_eps_tol_out*epsilon_check_a  && ...
    %         t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr) || ...
    %         (param.use_reweight_eps && t - reweight_last_step_iter > param.rw_tol)
    if (param.use_reweight_steps && t == reweight_steps(rw_counts) && t < param.reweight_max_reweight_itr) || ...
        (param.use_reweight_eps && rel_val(t) < param.reweight_rel_var && ...
        norm_residual_check <= param.adapt_eps_tol_out*norm_epsilon_check && ...
        t - reweight_last_step_iter > param.reweight_min_steps_rel_obj && t < param.reweight_max_reweight_itr) || ...
        (t - reweight_last_step_iter > 3000)
        
        fprintf('Reweighting: %i\n\n', reweight_step_count);
        for k = 1 : P
            d_val = abs(Psit{k}(xsol));
            weights1{k} = reweight_alpha ./ (reweight_alpha + d_val);
            %! [P.-A.] not used in the algo
            % weights1{k}(d_val > max(d_val) * param.reweight_abs_of_max) = 0;
        end

        % Calculate residual images
        res = zeros(size(xsol));
        for i = 1 : c
            Fx = A(xsol(:,:,i));
            g2 = zeros(No,1);
            for j = 1 : length(G{i})
                res_f = y{i}{j} - G{i}{j} * Fx(W{i}{j});
                u2{i}{j} = Gt{i}{j} * res_f;
                g2(W{i}{j}) = g2(W{i}{j}) + u2{i}{j};
            end
            res(:,:,i) = real(At(g2));
        end

        reweight_alpha = reweight_alpha_ff .* reweight_alpha;
        param.reweight_alpha = reweight_alpha;
        param.init_reweight_step_count = reweight_step_count+1;
        param.init_reweight_last_iter_step = t;
        param.init_t_start = t+1; 

        %! [P.-A.] "reweight_step_count > param.total_reweights" replaced by
        %! "reweight_step_count >= param.total_reweights"
        if (reweight_step_count >= param.total_reweights)
            param.reweight_max_reweight_itr = t+1;
            % fprintf('\n\n No more reweights \n\n');
        end
        
        if (reweight_step_count == 0) || (reweight_step_count == 1) || (~mod(reweight_step_count,5))
            % Save parameters (matfile solution)
            m = matfile([name, '_', ...
              num2str(param.ind) '_ch_' num2str(ch) '_'  num2str(param.gamma) '_rw=' num2str(reweight_step_count) '.mat'], ...
              'Writable', true);
            m.param = param;
            m.res = res;
            m.g = g;
            m.xsol = xsol;
            m.epsilon = epsilon;
            m.v2 = v2;
            m.proj = proj;
            m.t_block = t_block;
            m.norm_res = norm_res;
            m.v1 = v1;
            m.weights1 = weights1;
            m.res = res;
            % m.SNR = SNR;
            % m.SNR_average = SNR_average;
            m.l11 = l11;
            m.end_iter = end_iter;
            m.t_l11 = t_l11;
            m.t_master = t_master;
            m.t_data = t_data;
            m.rel_val = rel_val;
            clear m
            
            % save(['./results/result_SARA_cube_' num2str(param.ind) '_ch_' num2str(ch) '_gamma_' num2str(param.gamma) '_' num2str(reweight_step_count) '.mat'],'-v7.3','xsol','res', 'param');
            % Log
            if (param.verbose >= 1)
                fprintf('Backup iter: %i\n',t);
                fprintf('l11-norm = %e, rel_val = %e\n', l11(t), rel_val(t));
                fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
                fprintf(' SNR = %e\n', SNR);
            end
        end
       
        if (reweight_step_count >= param.total_reweights)
            fprintf('\n\n No more reweights \n\n');
            break;
            %! [P.-A.] why the instruction below? 
            % param.max_iter = t + 1000; 
        end

        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;
        rw_counts = rw_counts + 1;
    end
end

% Calculate residual images
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

m = matfile([name, '_', ...
    num2str(param.ind) '_ch_' num2str(ch) '_'  num2str(param.gamma) '_rw=' num2str(reweight_step_count) '.mat'], ...
'Writable', true);
m.param = param;
m.res = res;
m.g = g;
m.xsol = xsol;
m.epsilon = epsilon;
m.v2 = v2;
m.proj = proj;
m.t_block = t_block;
m.norm_res = norm_res;
m.v1 = v1;
m.weights1 = weights1;
m.res = res;

% Update param structure and save
param.reweight_alpha = reweight_alpha;
param.init_reweight_step_count = reweight_step_count;
param.init_reweight_last_iter_step = t;
param.init_t_start = t+1;
m.param = param;

% compute SNR
% m.SNR = SNR;
% m.SNR_average = SNR_average;
m.end_iter = end_iter;
m.t_l11 = t_l11;
m.t_master = t_master;
m.t_data = t_data;
m.rel_val = rel_val;
clear m

% Final log
if (param.verbose > 0)
    if (flag == 1)
        fprintf('Solution found\n');
        fprintf('Iter %i\n',t);
        fprintf(' Relative variation = %e\n', rel_val(t));
        fprintf(' Final residual = %e\n', residual_check);
        fprintf(' epsilon = %e\n', epsilon_check);
    else
        fprintf('Maximum number of iterations reached\n');
        fprintf('Iter %i\n',t);
        fprintf(' Relative variation = %e\n', rel_val(t));
        fprintf(' Final residual = %e\n', residual_check);
        fprintf(' epsilon = %e\n', epsilon_check);
    end
end
end

function [v1_, u1_, l11_, t_l11_] = run_par_waverec(v1_, Psit, Psi, xhat, weights1_, beta1)

tw = tic;
r1 = v1_ + Psit(xhat);
v1_ = r1 - sign(r1).*max(abs(r1) - beta1*weights1_,0);
u1_ = Psi(v1_);
t_l11_ = toc(tw);

% local L11 norm of current solution
l11_ = norm(r1(:),1);

end

function p = proj_l2ball(x, eps, y)
% projection of x onto the l2 ball centered in y with radius eps
p = x-y ;
p = p* min(eps/norm(p(:)),1) ;
p = p+y ;

end
