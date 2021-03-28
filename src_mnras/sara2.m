function [xsol,param,v1,v2,g,weights1,proj,t_block,reweighting_alpha,epsilon,t,rel_val,l11,norm_res,res,t_l11,t_master,end_iter] = ...
    sara2(y, epsilon, A, At, pU, G, W, Psi, Psit, param, init_file_name, name, x0, flag_homotopy, varargin)

% This function solves:
%
% min lambda * ||Psit(X)||_2,1   s.t.  || Y - A(X) ||_2 <= epsilon and x>=0
%
% Author: Abdullah Abdulaziz
c = size(y,2);
P = length(Psit);

% oversampling vectorized data length
No = size(W{1}{1}, 1);

% number of pixels
[M, N] = size(At(zeros(No, 1)));

%! ask Abdullah here
%maxNumCompThreads(12);
% util_create_pool(15);
% total number of workers (Q: facets workers, K: data workers)
numworkers = 12;
cirrus_cluster = parcluster('local');
cirrus_cluster.NumWorkers = numworkers;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
if cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads > ncores
    exit(1);
end
% maxNumCompThreads(param.num_workers);
parpool(cirrus_cluster, numworkers);
spmd
    dwtmode('zpd')
end

% Initializations
init_flag = isfile(init_file_name);
if init_flag
    init_m = matfile(init_file_name);
    fprintf('Resume from file %s\n\n', init_file_name)
end

%! -- TO BE CHECKED (primal initialization)
if init_flag
    xsol = init_m.xsol;
    param = init_m.param;
    epsilon = init_m.epsilon;
    fprintf('xsol, param and epsilon uploaded \n\n')
else
    if ~isempty(varargin)
        xsol = varargin{0};
    else
        xsol = zeros(M,N,c);
    end
    fprintf('xsol initialized \n\n')
end
xlast_reweight = xsol; %! assumes backup file exactly saved at the time a reweighting step occured
%! --

if init_flag
    g = init_m.g;
    fprintf('g uploaded \n\n')
else
    g = zeros(size(xsol));
    fprintf('g initialized \n\n')
end 

%! -- TO BE CHECKED
% Reweighting parameters
sig_bar = param.reweighting_sig_bar;
sig = param.reweighting_sig;
reweighting_alpha = param.reweighting_alpha;
reweighting_alpha_ff = param.reweighting_alpha_ff;

if isfield(param,'init_reweight_step_count')
    reweight_step_count = param.init_reweight_step_count;
    fprintf('reweight_step_count uploaded\n\n')
else
    reweight_step_count = 1;
    fprintf('reweight_step_count initialized \n\n')
end

if isfield(param,'init_reweight_last_iter_step')
    reweight_last_step_iter = param.init_reweight_last_iter_step;
    fprintf('reweight_last_iter_step uploaded \n\n')
else
    reweight_last_step_iter = 1;
    fprintf('reweight_last_iter_step initialized \n\n')
end

%! -- TO BE CHECKED (initialization of the weights)
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
        u1{k} = zeros(size(Psi{k}(v1{k})));
        % initial L1 descent step
        d_val = abs(Psit{k}(xsol));
        weights1{k} = reweighting_alpha ./ (reweighting_alpha + d_val);
        % weights1{k} = ones(size(v1{k},1),1);
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
    fprintf('t_start, t_block uploaded \n\n')
else
    for i = 1 : c
        t_block{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            t_block{i}{j} = 0;
        end
    end
    t_start = 1;
    fprintf('t_start, t_block initialized \n\n')
end

count_eps_update_down = 0;
count_eps_update_up = 0;

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
    l11 = NaN;
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
        fprintf('l11-norm = %e, rel_val = %e\n', l11, rel_val(t_start-1));
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

%%

fprintf('START THE LOOP MNRAS ver \n\n')

for t = t_start : param.reweighting_max_iter*param.pdfb_max_iter
    
    start_iter = tic;
    
    % update primal variable
    tw = tic;
    prev_xsol = xsol;
    xsol = max(real(xsol - g), 0);
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

            %! [P.-A.] what's the difference between residual_check_c and residual_check_a?
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
    l11 = sum(cell2mat(l11_cell));
    fprintf('Iter = %i, Time = %e, t_l11 = %e, t_data = %e\n',t,end_iter(t),t_data(t),t_l11(t));

    % Free memory
    Ftx=[]; g1=[];
    
    %% Display
    if ~mod(t,100)

        % SNR
        SNR = 20*log10(norm(x0(:))/norm(x0(:)-xsol(:)));
        psnrh = zeros(c,1);
        for i = 1:c
            psnrh(i) = 20*log10(norm(x0(:,i))/norm(x0(:,i)-reshape(xsol(:,:,i), [M*N, 1])));
        end
        SNR_average = mean(psnrh);
        
        % Log
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf('l11-norm = %e, rel_val = %e\n', l11, rel_val(t));
            %! [P.-A.] not sure why we have residual_c / residual_a, ...
            % fprintf(' epsilon = %e, residual = %e\n', norm(epsilon_check), norm(residual_check));
            % fprintf(' epsilon_c = %e, residual_c = %e\n', epsilon_check_c, residual_check_c);
            % fprintf(' epsilon_a = %e, residual_a = %e\n', epsilon_check_a, residual_check_a);
            fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
            fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
        end
    end

    %% Check convergence pdfb (inner solver)
    %! -- TO BE CHECKED
    pdfb_converged = (t - reweight_last_step_iter + 1 > param.pdfb_min_iter) && ...                                               % minimum number of pdfb iterations
        ( t - reweight_last_step_iter + 1 >= param.pdfb_max_iter || ...                                                          % maximum number of pdfb iterations reached
            (rel_val(t) < param.pdfb_rel_var && norm_residual_check <= param.pdfb_fidelity_tolerance*norm_epsilon_check) ... % relative variation and data fidelity within tolerance
        );
    
    %% Update epsilons
    flag_epsilonUpdate = param.use_adapt_eps && ...  % activate espilon update 
    (t > param.adapt_eps_start) && ...               % update allowed after a minimum of iterations in the 1st reweighting
    (rel_val(t) < param.adapt_eps_rel_var);          % relative variation between 2 consecutive pdfb iterations

    if flag_epsilonUpdate
        [epsilon, t_block] = update_epsilon(epsilon, t, t_block, ...
            norm_res, param.adapt_eps_tol_in, param.adapt_eps_tol_out, param.adapt_eps_steps, ...
            param.adapt_eps_change_percentage);
    end
    
    %% Reweighting
    if pdfb_converged
        rel_x_reweighting = norm(xlast_reweight(:) - xsol(:))/norm(xlast_reweight(:));
        xlast_reweight = xsol;

        fprintf('Reweighting: %i, relative variation: %e \n\n', reweight_step_count, rel_x_reweighting);

        reweighting_converged = pdfb_converged && ...                 % do not exit solver before the current pdfb algorithm converged
            reweight_step_count >  param.reweighting_min_iter && ...   % minimum number of reweighting iterations
            ( reweight_step_count >= param.reweighting_max_iter || ... % maximum number of reweighting iterations reached  
            rel_x_reweighting <= param.reweighting_rel_var ...        % relative variation
            );

        if reweighting_converged
            flag_convergence = 1;
            break;
        end

        % compute residual image
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

        % update weights
        for k = 1 : P
            d_val = abs(Psit{k}(xsol));
            weights1{k} = reweighting_alpha ./ (reweighting_alpha + d_val);
        end
        %! -- TO BE CHECKED
        if flag_homotopy
            reweighting_alpha = max(param.reweighting_alpha_ff * reweighting_alpha, 1);
        end
        %! --
        param.reweighting_alpha = reweighting_alpha;
        param.init_reweight_step_count = reweight_step_count+1;
        param.init_reweight_last_iter_step = t;
        param.init_t_start = t+1; 
        
        if (reweight_step_count == 1) || (reweight_step_count == 2) || (~mod(reweight_step_count,6))

            % compute SNR
            SNR = 20*log10(norm(x0(:))/norm(x0(:)-xsol(:)));
            psnrh = zeros(c,1);
            for i = 1:c
                psnrh(i) = 20*log10(norm(x0(:,i))/norm(x0(:,i)-reshape(xsol(:,:,i), [M*N, 1])));
            end
            SNR_average = mean(psnrh);

            % Save parameters (matfile solution)
            m = matfile([name, '_', ...
                num2str(param.cube_id), '_gam=',  num2str(param.gamma), '_rw=', num2str(reweight_step_count), '.mat'], ...
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

            % Log
            if (param.verbose >= 1)
                fprintf('Backup iter: %i\n',t);
                fprintf('l11-norm = %e, rel_val = %e\n', l11, rel_val(t));
                fprintf(' epsilon = %e, residual = %e\n', norm_epsilon_check, norm_residual_check);
                fprintf(' SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
            end
        end

        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;
        if (reweight_step_count >= param.reweighting_max_iter)
            fprintf('\n\n No more reweights \n\n');
        end
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
num2str(param.cube_id), '_gam=',  num2str(param.gamma), '_rw=', num2str(reweight_step_count), '.mat'], ...
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
param.reweighting_alpha = reweighting_alpha;
param.init_reweight_step_count = reweight_step_count;
param.init_reweight_last_iter_step = t;
param.init_t_start = t+1;
m.param = param;

% compute SNR
SNR = 20*log10(norm(x0(:))/norm(x0(:)-xsol(:)));
psnrh = zeros(c,1);
for i = 1:c
    psnrh(i) = 20*log10(norm(x0(:,i))/norm(x0(:,i)-reshape(xsol(:,:,i), [M*N, 1])));
end
SNR_average = mean(psnrh);
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
    if (flag_convergence == 1)
        fprintf('Solution found\n');
    else
        fprintf('Maximum number of iterations reached\n');
    end
    fprintf('Iter %i\n',t);
    fprintf(' L11-norm = %e, relative variation = %e\n', l11, rel_val(t));
    fprintf(' Final residual = %e\n', residual_check);
    fprintf(' epsilon = %e\n', epsilon_check);
    fprintf('SNR = %e, aSNR = %e\n\n', SNR, SNR_average);
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
