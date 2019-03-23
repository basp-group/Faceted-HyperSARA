function [xsol, L1_v, L1_vp, L2_v, L2_vp, delta_v, sol_v, snr_v, no_sub_itr_v, v1, v2, sol_reweight_v] = pdfb_bpcon_par_sing_sim_rescaled_precond(y, imsize, epsilont, epsilonts, epsilon, epsilons, A, At, T, pU, W, Psi, Psit, Psiw, Psitw, param)
%
% [xsol, L1_v, L1_vp, L2_v, L2_vp, delta_v, sol_v, snr_v, no_sub_itr_v, v1, v2, sol_reweight_v] = pdfb_bpcon_par_sim_rescaled_precond(y, epsilont, epsilonts, epsilon, epsilons, A, At, T, W, Psi, Psit, Psiw, Psitw, param)
% implements Algorithm 3 without randomisation as described in [1].
% The implementation simulates a distributed setup for parallel processing.
%
% Inputs:
% y{:} - the visibility data
% epsilont{:} - individual L2 bounds for each data split
% epsilonts{:} - individual L2 stop criterion bound for each data split
% epsilon - global L2 bound
% epsilons - global L2 stop criterion bound
% A - function handle, linear operator modeling the measurement process,
%     FFT and zero padding 
% At - the adjoint of A
% T{:} - gridding matrix containing the interpolation kernels (and any other
%     modeled DDEs)
% W{:} - selection of the FFT coefficients associated with each block of
%        data
% Psi - function handle, prior regularisation function
% Psit - the adjoint of Psi
% param - configuration parameters
% 
% Outputs: 
% xsol - final solution
% L1_v - evolution of the L1 norm
% L1_vp - evolution of the L1 norm per basis
% L2_v - evolution of the L2 norm
% L2_vp - evolution of the L2 norm per block
% delta_v - evolution of the relative solution variation
% sol_v - solutions at inerations provided in param.sol_steps
% snr_v - evolution of the SNR
% no_sub_itr_v - cell array containing cells with the number of sub-inteartions
%                performed at each iterations
% v1 - dual variable for the sparsity prior
% v2 - dual varaible for the data fidelity
% sol_reweight_v - solutions before performing reweighting steps
%
% Authors: Alexandru Onose
%
% [1] A. Onose, A. Dabbech, Y. Wiaux - An accelerated splitting algorithm 
% for radio-interferometric imaging: when natural and uniform weighting meet
% arXiv:1701.01748



% number of nodes
R = length(T);
P = length(Psit);

% number of pixels
Ny = imsize(1);
Nx = imsize(2);

%% optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end
if ~isfield(param, 'nu1')
    param.nu1 = ones(P, 1);
else
    if numel(param.nu1) == 1
        param.nu1 = ones(P, 1) * param.nu1;
    end
end
% if ~isfield(param, 'nu2')
%     param.nu2 = zeros(R, 1);
%     % maximum eigenvalue of operato A^T A
%     for q = 1:R
%         Tw = spalloc(size(T{q}, 1), No, size(T{q}, 2) * 16);
%         Tw(:, W{q}) = T{q};
%         fprintf('\nComputing operator norm: block %i \n', q)
%         param.nu2(q) = op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);
%         clear Tw;
%     end
%     fprintf('\n');
% else
if numel(param.nu2) == 1
    param.nu2 = ones(R, 1) * param.nu2;
end
% end
if ~isfield(param, 'sigma1'), param.sigma1 = 1./param.nu1; end
if ~isfield(param, 'sigma2'), param.sigma2 = 1./param.nu2; end
if ~isfield(param, 'gamma'), param.gamma = 1e-3; end
if ~isfield(param, 'tau'), param.tau = 0.49; end
if ~isfield(param, 'weights')
    param.weights = cell(P, 1);
    for k = 1:P
        param.weights{k} = ones(Nx * Ny, 1);
    end
else
    if ~iscell(param.weights)
        weights = param.weights;
        param.weights = cell(P, 1);
        for k = 1:P
            param.weights{k} = weights;
        end
    end
end
if isfield(param, 'initsol')
    xsol = param.initsol;
else
    % start from zero solution
%     xsol = At(zeros(No, 1));
    xsol = zeros(Ny, Nx);
end
if isfield(param, 'initv1')
    norm1 = cell(P, 1);
    r1 = cell(P, 1);
    u1 = cell(P, 1);
    v1 = cell(P, 1);
    if iscell(param.initv1)
        % initial dual variables
        v1 = param.initv1;
    else
        for k = 1:P
            % initial L1 part variable
            v1{k} = param.initv1;
        end
    end
    for k = 1:P
        % initial L1 descent step
        u1{k} = zeros(size(Psi{k}(v1{k})));
    end
    vy1 = v1;
else
    norm1 = cell(P, 1);
    r1 = cell(P, 1);
    u1 = cell(P, 1);
    v1 = cell(P, 1);
    for k = 1:P
        % start from zero solution
        v1{k} = zeros(size(Psit{k}(xsol)));
        
        % initial L1 descent step
        u1{k} = zeros(size(Psi{k}(v1{k})));
    end
    vy1 = v1;
end
if isfield(param, 'initv2')
    norm2 = cell(R, 1);
    r2 = cell(R, 1);
    u2 = cell(R, 1);
    v2 = cell(R, 1);
    % initial L2 ball variables
    % initial dual variables
    if iscell(param.initv2)
        % initial dual variables
        v2 = param.initv2;
    else
        for q = 1:R
            % initial L1 part variable
            v2{q} = param.initv2;
        end
    end
    for q = 1:R
        % initial L1 part descent step
        u2{q} = zeros(size(T{q}, 1), 1);
    end
    vy2 = v2;
else
    norm2 = cell(R, 1);
    r2 = cell(R, 1);
    u2 = cell(R, 1);
    v2 = cell(R, 1);
    for q = 1:R
        % initial L1 part variable
        v2{q} = zeros(length(y{q}), 1);
        % initial L1 part descent step
        u2{q} = zeros(size(T{q}, 2), 1);
    end
    vy2 = v2;
end
if ~isfield(param, 'lambda0'), param.lambda0 = 1; end
if ~isfield(param, 'lambda1'), param.lambda1 = 1; end
if ~isfield(param, 'lambda2'), param.lambda2 = 1; end
if ~isfield(param, 'sol_steps')
    param.sol_steps = inf;
else
    if param.sol_steps(end) ~= inf
        param.sol_steps = [param.sol_steps inf];
    end
end
if ~isfield(param, 'reweight_steps')
    param.reweight_steps = inf;
else
    if param.reweight_steps(end) ~= inf
        param.reweight_steps = [param.reweight_steps inf];
    end
end
if ~isfield(param, 'best_bound_steps')
    param.best_bound_steps = inf;
else
    if param.best_bound_steps(end) ~= inf
        param.best_bound_steps = [param.best_bound_steps inf];
    end
end
if ~isfield(param, 'omega1'), param.omega1 = 1; end
if ~isfield(param, 'omega2'), param.omega2 = 1; end
if ~isfield(param, 'im0')
    param.im0 = zeros(Ny, Nx);
end
if ~isfield(param, 'global_stop_bound'), param.global_stop_bound = 1; end

if ~isfield(param, 'use_proj_elipse_fb')
    param.use_proj_elipse_fb = 1;
end
if ~isfield(param, 'elipse_proj_max_iter')
    param.elipse_proj_max_iter = 2000;
end
if ~isfield(param, 'elipse_proj_min_iter')
    param.elipse_proj_min_iter = 1;
end
if ~isfield(param, 'elipse_proj_eps')
    param.elipse_proj_eps = 1e-8;
end
if ~isfield(param, 'use_reweight_steps')
    param.use_reweight_steps = 0;
end
if ~isfield(param, 'use_reweight_eps')
    param.use_reweight_eps = 0;
end
if ~isfield(param, 'use_best_bound_steps')
    param.use_best_bound_steps = 0;
end
if ~isfield(param, 'use_best_bound_eps')
    param.use_best_bound_eps = 0;
end


%% set up log variables
L1_v = zeros(param.max_iter, 1);
L1_vp = zeros(param.max_iter, P);
L2_v = zeros(param.max_iter, 1);
L2_vp = zeros(param.max_iter, R);
no_sub_itr_v = cell(param.max_iter, 1);

delta_v = zeros(param.max_iter, 1);

sol_steps = param.sol_steps;
sol_step_count = 1;

reweight_steps = param.reweight_steps;
reweight_step_count = 1;
reweight_last_step_iter = 1;

sol_v = zeros(length(sol_steps)-1, Ny, Nx);
sol_reweight_v = zeros(0, Ny, Nx);

snr_v = zeros(param.max_iter, 1);

%% useful functions for the projection
% scaling, projection on L2 norm
sc = @(z, radius) z * min(radius/norm(z(:)), 1);

% thresholding negative values
hardt = @(z) max(real(z), min(-param.im0, 0));
% hardt = @(z) max(real(z), 0);

%soft thresholding operator
soft = @(z, T) sign(z) .* max(abs(z)-T, 0); 


%% initialization

% initial primal gradient like step
g1 = zeros(size(xsol));
g2 = zeros(size(xsol));

% solution flag: 0 - max iteration reached; 1 - solution found
flag = 0;

%% store useful variables
% step size for the dual variables
sigma1 = param.sigma1;
sigma2 = param.sigma2;

% step size primal 
tau = param.tau;

% relaxation parameters
lambda0 = param.lambda0;
lambda1 = param.lambda1;
lambda2 = param.lambda2;


% weights
weights = param.weights;

% omega sizes
omega1 = param.omega1;
omega2 = param.omega2;

gamma = param.gamma;

reweight_alpha = param.reweight_alpha;
reweight_alpha_ff = param.reweight_alpha_ff;

% [proj] = solver_find_elipse_point(y, pU, A, T, xsol, v2, W, epsilont, param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);


%% main loop: sequential + simulated parallel
% xsol     - current solution estimate
% prev_sol - previous solution estimate
% r1       - L1 part of criterion: r1 = Psi' * sol
% v1       - L1 dual variable
% r2       - L2 part of criterion: r2 = T * A * sol
% v2       - L2 dual variable

for t = 1:param.max_iter    
    tm = tic;
    
    %% primal update
    ysol = hardt(xsol - tau*(omega1 * g1 + omega2 * g2));
    prev_xsol = xsol;
    xsol = xsol + lambda0 * (ysol - xsol);
    
    norm_prevsol = norm(prev_xsol(:));
    % solution relative change
    if (norm_prevsol == 0)
        rel_sol_norm_change = 1;
    else
        rel_sol_norm_change = norm(xsol(:) - prev_xsol(:))/norm_prevsol;
    end
    
    prev_xsol = 2*ysol - prev_xsol;
    
    %% L1 prox update: dual variables update
    % parallel for all bases
    for k = 1:P        
        r1{k} = weights{k} .* Psit{k}(prev_xsol);
        vy1{k} = v1{k} + r1{k} - soft(v1{k} + r1{k}, gamma / sigma1(k));
        
        v1{k} = v1{k} + lambda1 * (vy1{k} - v1{k});
        u1{k} = Psi{k}(weights{k} .* v1{k});
        
        % local L1 norm of current solution
        norm1{k} = sum(abs(r1{k}));
    end

    
    %% L2 ball projection update: dual variables update
    
    % non gridded measurements of curent solution 
    ns = A(prev_xsol);
    
    % partial non gridded measurements for each node
    ns_p = cell(R, 1);
   
    % select parts to be sent to nodes
    for q = 1:R
        ns_p{q} = ns(W{q});
    end
    
    no_sub_itr = cell(R, 1);
    res = cell(R, 1);
    
    % parallel for all R blocks
    for q = 1:R
        r2{q} = T{q} .* ns_p{q};
        
        if param.use_proj_elipse_fb
            fprintf('Preconditioning step');
            if t == 1
                proj{q} = zeros(size(y{q}));
            end
            [proj{q}, no_sub_itr{q}] = solver_proj_elipse_fb(1 ./ pU{q} .* v2{q}, r2{q}, y{q}, pU{q}, epsilont{q}, proj{q}, param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);
            vy2{q} = v2{q} + pU{q} .* r2{q} - pU{q} .* proj{q};
        else
            vy2{q} = v2{q} + r2{q} - y{q} - sc(v2{q} + r2{q} - y{q}, epsilont{q});
        end
%         if ~param.use_proj_elipse_fb
%             [proj{q}, no_sub_itr{q}] = solver_proj_elipse((sqrt(pU{q}) .* 1./pU{q}) .* (v2{q} + pU{q} .* r2{q}), 1 ./ sqrt(pU{q}), epsilont{q}, proj{q}, y{q}, param.elipse_proj_max_iter, param.elipse_proj_eps);
%             vy2{q} = v2{q} + pU{q} .* r2{q} - (pU{q} .* 1 ./ sqrt(pU{q})) .* proj{q};
%             [proj{q}, no_sub_itr{q}] = solver_proj_elipse(1 ./ sqrt(pU{q}) .* (v2{q} + pU{q} .* r2{q}), 1 ./ sqrt(pU{q}), epsilont{q}, proj{q}, y{q}, param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);
%             vy2{q} = v2{q} + pU{q} .* r2{q} - sqrt(pU{q}) .* proj{q};
%         end
          

        v2{q} = v2{q} + lambda2 * (vy2{q} - v2{q});
        u2{q} = T{q} .* v2{q};

        % norm of residual
        res{q} = r2{q} - y{q};
        norm2{q} = norm(res{q});
    end


    %% update the primal gradient
    g1 = zeros(size(xsol));
    g2 = zeros(size(xsol));
    
    for k = 1:P
        g1 = g1 + sigma1(k) * u1{k};
    end
    
    uu = zeros(Ny*Nx, 1);
    for q = 1:R
        uu(W{q}) = sigma2(q) * u2{q};
    end
    g2 = g2 + At(uu);
    
    tm = toc(tm);
          
    %% stopping criterion and logs
 
    % log
    if (param.verbose >= 1)
        fprintf('Iter %i\n',t);
        fprintf(' L1 norm                       = %e\n', sum(cell2mat(norm1)));
        fprintf(' Residual                      = %e\n', norm(cell2mat(norm2)));
        fprintf(' Global residual bound         = %e\n', epsilon);
        fprintf(' Distributed residual L2 ball  = %e\n', norm(cell2mat(epsilont)));
        fprintf(' Distributed residual L2 bound = %e\n', norm(cell2mat(epsilonts)));
        fprintf(' Relative solution norm change = %e\n\n', rel_sol_norm_change);
        
        if (param.verbose >= 2)
            for q = 1:R
                fprintf('   Residual %i                     = %e\n', q, norm2{q});
                fprintf('   Residual L2 ball %i             = %e\n', q, epsilont{q});
                fprintf('   Residual L2 bound %i            = %e\n\n', q, epsilonts{q});
            end
        end
        
        fprintf('Time for iteration %i: %3.3f\n',t, tm);
    end
    if (param.verbose <= 0.5)
        fprintf('.\n');fprintf('\b');
        if mod(t, 50) == 0
            fprintf('\n');
        end
    end
    if (param.verbose >= 0.5)
        L1_v(t) = sum(cell2mat(norm1));
        L2_v(t) = norm(cell2mat(norm2));
        no_sub_itr_v{t} = no_sub_itr;
        for q = 1:R
            L2_vp(t, q) = norm2{q};
        end
        for k = 1:P
            L1_vp(t, k) = norm1{k};
        end
        delta_v(t) = rel_sol_norm_change;
        try 
            snr_v(t) = 20*log10(norm(param.im(:))/norm(param.im(:) - xsol(:)));
        catch
            snr_v(t) = 0;
        end
        
        if t == sol_steps(sol_step_count)
            sol_v(sol_step_count, :, :) = xsol;
            sol_step_count = sol_step_count + 1;
        end
    end
    
    
    if (param.use_reweight_steps && t == reweight_steps(reweight_step_count)) || ...
        (param.use_reweight_eps && ...
        norm(cell2mat(norm2)) <= norm(cell2mat(epsilonts)) && ...
        param.reweight_min_steps_rel_obj < t - reweight_last_step_iter && ...
        rel_sol_norm_change < param.reweight_rel_obj)
    
        % parallel for all bases
        for k = 1:P
            d_val = abs(Psit{k}(xsol));
            weights{k} = reweight_alpha ./ (reweight_alpha + d_val);
        end

        reweight_alpha = reweight_alpha_ff * reweight_alpha;
        weights_mat = [weights{:}];
        weights_mat = weights_mat(:);
        sigma1(:) = 1/op_norm(@(x) weights_mat .* Psitw(x), @(x) Psiw(weights_mat .* x), [Ny, Nx], 1e-8, 200, 0);
        
        fprintf('\n\n\n\n\n\n\n Performed reweight no %d \n\n\n\n\n', reweight_step_count);
        
        sol_reweight_v(reweight_step_count, :, :) = xsol;
        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;
    end

    % global stopping criteria
    if rel_sol_norm_change < param.rel_obj && ...
            ((param.global_stop_bound && norm(cell2mat(norm2)) <= norm(cell2mat(epsilonts))) || ...
            (~param.global_stop_bound && prod(cell2mat(norm2) <= cell2mat(epsilonts))))
        flag = 1;
        break;
    end
end


% final log
if (param.verbose > 0)
    if (flag == 1)
        fprintf('\nSolution found\n');
        fprintf(' L1 norm                       = %e\n', sum(cell2mat(norm1)));
        fprintf(' Residual                      = %e\n', norm(cell2mat(norm2)));
        fprintf(' Global residual bound         = %e\n', epsilon);
        fprintf(' Distributed residual L2 ball  = %e\n', norm(cell2mat(epsilont)));
        fprintf(' Distributed residual L2 bound = %e\n', norm(cell2mat(epsilonts)));
        fprintf(' Relative solution norm change = %e\n\n', rel_sol_norm_change);
        
        for q = 1:R
            fprintf('   Residual %i                     = %e\n', q, norm2{q});
            fprintf('   Residual L2 ball %i             = %e\n', q, epsilont{q});
            fprintf('   Residual L2 bound %i            = %e\n\n', q, epsilonts{q});
        end
    else
        fprintf('\nMaximum number of iterations reached\n');
        fprintf(' L1 norm                       = %e\n', sum(cell2mat(norm1)));
        fprintf(' Residual                      = %e\n', norm(cell2mat(norm2)));
        fprintf(' Global residual bound         = %e\n', epsilon);
        fprintf(' Distributed residual L2 ball  = %e\n', norm(cell2mat(epsilont)));
        fprintf(' Distributed residual L2 bound = %e\n', norm(cell2mat(epsilonts)));
        fprintf(' Relative solution norm change = %e\n\n', rel_sol_norm_change);
        
        for q = 1:R
            fprintf('   Residual %i                     = %e\n', q, norm2{q});
            fprintf('   Residual L2 ball %i             = %e\n', q, epsilont{q});
            fprintf('   Residual L2 bound %i            = %e\n\n', q, epsilonts{q});
        end
    end
end
fprintf('\n');

xsol = hardt(xsol);

sol_reweight_v(reweight_step_count, :, :) = xsol;        


% trim the log vectors to the size of the actual iterations performed
if (param.verbose >= 0.5)
    L1_v = L1_v(1:t);
    L1_vp = L1_vp(1:t, :);
    L2_v = L2_v(1:t);
    L2_vp = L2_vp(1:t, :);
    delta_v = delta_v(1:t);
    sol_v = sol_v(1:sol_step_count-1, :, :);
    snr_v = snr_v(1:t);
    no_sub_itr_v = {no_sub_itr_v{1:t}}';
end

end

