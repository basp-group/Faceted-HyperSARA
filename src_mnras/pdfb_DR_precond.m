function [xsol, t_block, epsilon, t, rel_fval, norm2, res, end_iter] = pdfb_DR_precond(y, epsilon, A, At, H, W, pU, T, Wm, Psi, Psit, param, reduction_version, realdatablocks, fouRed_gamma, typeStr, wterm, levelG, levelC, init_file_name)
%xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,res,end_iter
% Inputs:
% y{:} - the visibility data
% imsize - image size
% epsilon - L2 bound
% A - function handle, linear operator modeling the measurement process,
%     FFT and zero padding 
% At - the adjoint of A
% H{:} - holographic matrix
% pU{:} - pre-conditioning matrix
% T{:} - selected "singular values" for Fourier reduction
% W{:} - selection matrix for Fourier reduction
% Psi - function handle, prior regularisation function
% Psit - the adjoint of Psi
% param - configuration parameters
% reduction_version: 1. embedding operator F*Phi^t 2. embedding operator G^t 
% realdatablocks (for real data only): 2. 2 blocks 9. 9 blocks
% 
% Outputs: 
% xsol - final solution
% rel_fval - evolution of the relative solution variation
%
% Author: Ming Jiang, E-mail: ming.jiang@epfl.ch
%

FT2 = @(x) fftshift(fft2(ifftshift(x))) / sqrt(numel(x));
IFT2 = @(x) fftshift(ifft2(ifftshift(x))) * sqrt(numel(x));

monitor_step = 500;
reweight_monitor_step = 5;

% number of nodes
R = length(H);
P = length(Psit);

% Variable flag for the case where W is present
flagW = 0;
if ~isempty(W)
    flagW = 1;
end

% number of over-sampled pixels
if flagW
    No = size(W{1}, 1);
else
    No = size(H{1}, 2);
end

% number of pixels (spatial dimensions)
[Ny, Nx] = size(At(zeros(No, 1)));

% fprintf('rw_tol: %i \n\n', param.rw_tol)
% Initializations.
init_flag = isfile(init_file_name);
if init_flag
    init_m = matfile(init_file_name);
end

if init_flag
    fprintf('Load data from %s\n', init_file_name)
    param = init_m.param;
    xsol = init_m.xsol;
    epsilon = init_m.epsilon;
    weights = init_m.weights;
    norm1 = init_m.norm1;
    u1 = init_m.u1;
    v1 = init_m.v1;
    norm2 = init_m.norm2;
    r2 = init_m.r2;
    u2 = init_m.u2;
    v2 = init_m.v2;
    vy2 = v2;
    proj = init_m.proj;
    t_block = init_m.t_block;
    reweight_last_step_iter = param.init_reweight_last_iter_step;
    reweight_step_count = param.init_reweight_step_count;
    count_eps_update_down = init_m.count_eps_update_down;
    count_eps_update_up = init_m.count_eps_update_up;
    rel_fval = init_m.rel_fval;
    end_iter = init_m.end_iter;
    g1 = init_m.g1;
    g2 = init_m.g2;
    flag = init_m.flag;
    t_start = param.init_t_start;
    reweight_alpha= init_m.reweight_alpha;
%     peak_psf = init_m.peak_psf;
    t_master = init_m.t_master;
    t_l11 = init_m.t_l11;
    t_data = init_m.t_data;
else
    %% optional input arguments
    % if ~isfield(param, 'nu1')
    %     param.nu1 = ones(P, 1);
    % else
    %     if numel(param.nu1) == 1
    %         param.nu1 = ones(P, 1) * param.nu1;
    %     end
    % end

    % if numel(param.nu2) == 1
    %     param.nu2 = ones(R, 1) * param.nu2;
    % end
    % end
    if ~isfield(param, 'sigma1'), param.sigma1 = 1./param.nu1; end
    if ~isfield(param, 'sigma2'), param.sigma2 = 1./param.nu2; end
    if ~isfield(param, 'gamma'), param.gamma = 1e-3; end
    if ~isfield(param, 'tau'), param.tau = 0.49; end
    if isfield(param,'init_weights')
        weights = param.init_weights;
        fprintf('weights uploaded \n\n')
    else
        weights = cell(P, 1);
        for k = 1:P
            weights{k} = ones(Nx * Ny,1);
        end
        fprintf('weights NOT uploaded \n\n')
    end
    if isfield(param, 'initsol')
        xsol = param.initsol;
        fprintf('xsol uploaded \n\n')
    else
        % start from zero solution
        xsol = zeros(Ny, Nx);
        fprintf('xsol NOT uploaded \n\n')
    end
    if isfield(param, 'initv1')
        norm1 = cell(P, 1);
    %     r1 = cell(P, 1);
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
    %     vy1 = v1;
        fprintf('v1 uploaded \n\n')
    else
        norm1 = cell(P, 1);
    %     r1 = cell(P, 1);
        u1 = cell(P, 1);
        v1 = cell(P, 1);
        for k = 1:P
            % start from zero solution
            v1{k} = zeros(size(Psit{k}(xsol)));

            % initial L1 descent step
            u1{k} = zeros(size(Psi{k}(v1{k})));
        end
    %     vy1 = v1;
        fprintf('v1 NOT uploaded \n\n')
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
        fprintf('v2 uploaded \n\n')
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
        fprintf('v2 NOT uploaded \n\n')
    end

    if isfield(param,'init_t_block')
        t_block = param.init_t_block;
        reweight_last_step_iter = param.init_t;
        reweight_step_count = param.reweight_step_count+1;
        fprintf('t t_block uploaded \n\n')
    else
        t_block = cell(R, 1);
        for q = 1 : R
            t_block{q} = 0;
        end
        reweight_last_step_iter = 1;
        reweight_step_count = 0;
        fprintf('t t_block NOT uploaded \n\n')
    end
    %% set up log variables
    count_eps_update_down = 0;
    count_eps_update_up = 0;

    rel_fval = zeros(param.max_iter, 1);
    end_iter = zeros(param.max_iter, 1);
    t_master = zeros(param.max_iter, 1);
    t_l11 = zeros(param.max_iter, 1);
    t_data = zeros(param.max_iter, 1);
    
%     %% initialization
%     % psf
%     dirac=zeros(Ny,Nx);
%     dirac(Ny/2,Nx/2)=1;
%     dirac_uv=A(dirac);
%     for q = 1:R        
%         if reduction_version == 1
%             if flagW
%                 tmp = FT2(real(At(H{q} * dirac_uv(W{q}))));
%             else
%                 tmp = FT2(real(At(H{q} * dirac_uv)));
%             end
%             tmp = tmp(:);
%             tmp2{q} = T{q} .* tmp(Wm{q});
%         elseif reduction_version == 2
%             if flagW
%                 tmp2{q} = T{q} .* (H{q} * dirac_uv(W{q}));
%             else
%                 tmp2{q} = T{q} .* (H{q} * dirac_uv);
%             end
%         end
%     end
%     
%     g_f2 = zeros(No, 1);
%     for q = 1 : R
%         if reduction_version == 1   % H is self-adjoint in this case
%             tmp = zeros(size(Wm{q}));
%             tmp(Wm{q}) = T{q} .* tmp2{q};
%             if flagW
%                 g_f2(W{q}) = g_f2(W{q}) + H{q} * A(real(IFT2(reshape(tmp, Ny, Nx))));
%             else
%                 g_f2 = g_f2 + H{q} * A(real(IFT2(reshape(tmp, Ny, Nx))));
%             end
%         elseif reduction_version == 2
%             if flagW
%                 g_f2(W{q}) = g_f2(W{q}) + H{q}' * (T{q} .* tmp2{q});
%             else
%                 g_f2 = g_f2 + H{q}' * (T{q} .* tmp2{q});
%             end
%         end
%     end
%     psf = real(At(g_f2));
%     peak_psf = max(psf(:));
%     fprintf('Peak psf: %e\n', peak_psf)

    % initial primal gradient like step
    g1 = zeros(size(xsol));
    g2 = zeros(size(xsol));
    
    reweight_alpha = param.reweight_alpha;

    % solution flag: 0 - max iteration reached; 1 - solution found
    flag = 0;
    t_start = 1;
end

%% useful functions for the projection
% scaling, projection on L2 norm
sc = @(z, radius) z * min(radius/norm(z(:)), 1);

% thresholding negative values
% hardt = @(z) max(real(z), min(-param.im0, 0));
hardt = @(z) max(real(z), 0);

%soft thresholding operator
soft = @(z, T) sign(z) .* max(abs(z)-T, 0); 

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
% weights = param.weights;
reweight_steps = param.reweight_steps;

% omega sizes
omega1 = param.omega1;
omega2 = param.omega2;

gamma = param.gamma;

reweight_alpha_ff = param.reweight_alpha_ff;
reweight_abs_of_max = param.reweight_abs_of_max;

A = afclean(A);
At = afclean(At);
for k = 1 : P
    Psi{k} = afclean(Psi{k});
    Psit{k} = afclean(Psit{k});
end

%% main loop: sequential + simulated parallel
% xsol     - current solution estimate
% prev_sol - previous solution estimate

util_create_pool(P+1);

start_loop = tic;

for t = t_start:param.max_iter    
    tm = tic;
    
    %% primal update
    tw = tic;
    ysol = hardt(xsol - tau*(omega1 * g1 + omega2 * g2));
    prev_xsol = xsol;
    xsol = xsol + lambda0 * (ysol - xsol);
    t_master(t) = toc(tw);
    
    norm_prevsol = norm(prev_xsol(:));
    % solution relative change
    if (norm_prevsol == 0)
        rel_sol_norm_change = 1;
    else
        rel_sol_norm_change = norm(xsol(:) - prev_xsol(:))/norm_prevsol;
    end
    
    rel_fval(t) = rel_sol_norm_change;
    prev_xsol = 2*ysol - prev_xsol;

    %% L1 prox update: dual variables update
    % parallel for all bases
%     parfor k = 1:P        
%         r1{k} = weights{k} .* Psit{k}(prev_xsol);
%         vy1{k} = v1{k} + r1{k} - soft(v1{k} + r1{k}, gamma / sigma1);
%         
%         v1{k} = v1{k} + lambda1 * (vy1{k} - v1{k});
%         u1{k} = Psi{k}(weights{k} .* v1{k});
%         
%         % local L1 norm of current solution
%         norm1{k} = sum(abs(r1{k}));
%     end
    for k = 1:P
        f(k) = parfeval(@run_par_waverec, 4, v1{k}, Psit{k}, Psi{k}, prev_xsol, gamma, weights{k}, sigma1, lambda1);
    end

    %% L2 ball projection update: dual variables update
    
    % non gridded measurements of curent solution 
    ns = A(prev_xsol);

    res = cell(R, 1);
    
    % parallel for all R blocks
    residual_check_c = 0;
    epsilon_check_c = 0;
    residual_check_a = 0;
    epsilon_check_a = 0;
    tw = tic;
    for q = 1:R        
        if reduction_version == 1
            if flagW
                tmp = FT2(real(At(H{q} * ns(W{q}))));
            else
                tmp = FT2(real(At(H{q} * ns)));
            end
            tmp = tmp(:);
            r2{q} = T{q} .* tmp(Wm{q});
        elseif reduction_version == 2
            if flagW
                r2{q} = T{q} .* (H{q} * ns(W{q}));
            else
                r2{q} = T{q} .* (H{q} * ns);
            end
        end
        
        if param.use_proj_elipse_fb
%             fprintf('Preconditioning step\n');
            if t == t_start
                proj{q} = zeros(size(y{q}));
            end
            [proj{q}, ~] = solver_proj_elipse_fb(1 ./ pU{q} .* v2{q}, r2{q}, y{q}, pU{q}, epsilon{q}, proj{q}, param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);
            vy2{q} = v2{q} + pU{q} .* r2{q} - pU{q} .* proj{q};
        else
            vy2{q} = v2{q} + r2{q} - y{q} - sc(v2{q} + r2{q} - y{q}, epsilon{q});
        end
        v2{q} = v2{q} + lambda2 * (vy2{q} - v2{q});
        u2{q} = T{q} .* v2{q};

        % norm of residual
        res{q} = y{q} - r2{q};
        norm2{q} = norm(res{q});
        
%         residual_check = norm2{q};
%         epsilon_check = epsilon{q};
        
        if realdatablocks == 2
            if q == 1
                residual_check_c = norm2{q};
                epsilon_check_c = epsilon{q};
            else
                residual_check_a = norm2{q};
                epsilon_check_a = epsilon{q};
            end
            
        elseif realdatablocks == 9
            if q == 1 || q == 2
                residual_check_c = sqrt(residual_check_c + norm2{q}^2);
                epsilon_check_c = sqrt(epsilon_check_c + epsilon{q}^2);
            else
                residual_check_a = sqrt(residual_check_a + norm2{q}^2);
                epsilon_check_a = sqrt(epsilon_check_a + epsilon{q}^2);
            end
        end      
    end
    t_data(t) = toc(tw);
    residual_check = sqrt(residual_check_a^2 + residual_check_c^2);
    epsilon_check = sqrt(epsilon_check_a^2 + epsilon_check_c^2);
    
%     if realdatablocks == 9
%         epsilon_check_c = sqrt(epsilon_check_c);
%         residual_check_c = sqrt(residual_check_c);
% 
%         epsilon_check_a = sqrt(epsilon_check_a);
%         residual_check_a = sqrt(residual_check_a);
%     end
    
    %% ADAPTIVE bound update on each block
    if( param.use_adapt_eps ==1) && t > param.adapt_eps_start
     for q = 1:R
       if (norm2{q} < param.adapt_eps_tol_in * epsilon{q}) && (t > t_block{q} + param.adapt_eps_steps) && ...
               (rel_sol_norm_change <  param.adapt_eps_rel_obj)
           t_block{q} = t;
           epsilon{q} = norm2{q} + (-norm2{q} + epsilon{q}) * (1 - param.adapt_eps_change_percentage);
           count_eps_update_down = count_eps_update_down +1;
           fprintf ('Updated  epsilon DOWN: %e\t, residual: %e\t, Block: %i\n', epsilon{q}, norm2{q}, q); 
       end

       if (norm2{q} > param.adapt_eps_tol_out * epsilon{q}) && (t > t_block{q} + param.adapt_eps_steps) && ...
               (rel_sol_norm_change < param.adapt_eps_rel_obj)
           t_block{q} = t;
%            target_eps = norm2{q} + (norm2{q} - epsilon{q}) * param.adapt_eps_change_percentage;
%            if target_eps > param.l2_upper_bound{q}
%                epsilon{q} = param.l2_upper_bound{q};
%            else
%                epsilon{q} = target_eps;
%            end
           epsilon{q} = norm2{q} + (norm2{q} - epsilon{q}) * param.adapt_eps_change_percentage;
           count_eps_update_up = count_eps_update_up +1;
           fprintf('Updated epsilon UP: %e\t, residual: %e\t, Block: %i\n', epsilon{q}, norm2{q}, q);
       end
     end
    end 


    %% update the primal gradient
%     g1 = sparse(Ny, Nx);
    g1 = zeros(size(xsol));
    
%     for k = 1:P
%         g1 = g1 + sigma1 * u1{k};
%     end
    
    for k = 1:P
        [idx, v1_, u1_, l1_, t_l11_] = fetchNext(f);
        v1{idx} = v1_;
        u1{idx} = u1_;
        norm1{idx} = l1_;
        t_l11(t) = t_l11(t) + t_l11_;
        g1 = g1 + u1{idx};
    end
    
    uu = zeros(No, 1);
    for q = 1:R
        if reduction_version == 1
            tmp = zeros(size(Wm{q}));
            tmp(Wm{q}) = u2{q};
            if flagW
                uu(W{q}) = uu(W{q}) + H{q} * A(real(IFT2(reshape(tmp, Ny, Nx))));
            else
                uu = uu + H{q} * A(real(IFT2(reshape(tmp, Ny, Nx))));
            end                
        elseif reduction_version == 2
            if flagW
                uu(W{q}) = uu(W{q}) + H{q}' * u2{q};
            else
                uu = uu + H{q}' * u2{q};
            end
        end
    end
    g2 = real(At(sigma2 * uu));   
    end_iter(t) = toc(tm);
    t_l11(t) = t_l11(t)/P; % average compute time for the dual variable
    fprintf('Iter = %i, Time = %e, t_master = %e, t_l11 = %e, t_data = %e, Rel_error = %e \n',t, end_iter(t), t_master(t), t_data(t), t_l11(t), rel_fval(t));      
    %% stopping criterion and logs
 
    % log
    if ~mod(t,monitor_step)
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf(' L1 norm                       = %e\n', sum(cell2mat(norm1)));
            fprintf(' Residual                      = %e\n', norm(cell2mat(norm2)));
            fprintf(' Epsilon                       = %e\n', norm(cell2mat(epsilon)));
            fprintf(' Relative solution norm change = %e\n\n', rel_sol_norm_change);

            if (param.verbose >= 2)
                for q = 1:R
                    fprintf('   Residual %i                     = %e\n', q, norm2{q});
                    fprintf('   Epsilon %i             = %e\n', q, epsilon{q});
                end
            end
        end
        
        % Calculate residual image
        g_f2 = zeros(No, 1);
        for q = 1 : R
            if reduction_version == 1   % H is self-adjoint in this case
                tmp = zeros(size(Wm{q}));
                tmp(Wm{q}) = T{q} .* res{q};
                if flagW
                    g_f2(W{q}) = g_f2(W{q}) + H{q} * A(real(IFT2(reshape(tmp, Ny, Nx))));
                else
                    g_f2 = g_f2 + H{q} * A(real(IFT2(reshape(tmp, Ny, Nx))));
                end
            elseif reduction_version == 2
                if flagW
                    g_f2(W{q}) = g_f2(W{q}) + H{q}' * (T{q} .* res{q});
                else
                    g_f2 = g_f2 + H{q}' * (T{q} .* res{q});
                end
            end
        end
        res_im = real(At(g_f2));
        
        if wterm
            prev_file = ['results/xsol_ch', num2str(param.chInd), '_prec_it', num2str(t-monitor_step), '_gamma', num2str(gamma), '_', num2str(realdatablocks),...
                'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '.fits'];
            if isfile(prev_file)
                delete(prev_file)
            end
            fitswrite(xsol, ['results/xsol_ch', num2str(param.chInd), '_prec_it', num2str(t), '_gamma', num2str(gamma), '_', num2str(realdatablocks),...
                'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '.fits']);
            
            prev_file = ['results/res_ch', num2str(param.chInd), '_prec_it', num2str(t-monitor_step), '_gamma', num2str(gamma), '_', num2str(realdatablocks),...
                'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '.fits'];
            if isfile(prev_file)
                delete(prev_file)
            end
            fitswrite(res_im, ['results/res_ch', num2str(param.chInd), '_prec_it', num2str(t), '_gamma', num2str(gamma), '_', num2str(realdatablocks),...
                'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '.fits']);
        else
            prev_file = ['results/xsol_ch', num2str(param.chInd), '_prec_it', num2str(t-monitor_step), '_gamma', num2str(gamma), '_', num2str(realdatablocks),...
                'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.fits'];
            if isfile(prev_file)
                delete(prev_file)
            end
            fitswrite(xsol, ['results/xsol_ch', num2str(param.chInd), '_prec_it', num2str(t), '_gamma', num2str(gamma), '_', num2str(realdatablocks),...
                'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.fits']);
            
            prev_file = ['results/res_ch', num2str(param.chInd), '_prec_it', num2str(t-monitor_step), '_gamma', num2str(gamma), '_', num2str(realdatablocks),...
                'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.fits'];
            if isfile(prev_file)
                delete(prev_file)
            end
            fitswrite(res_im, ['results/res_ch', num2str(param.chInd), '_prec_it', num2str(t), '_gamma', num2str(gamma), '_', num2str(realdatablocks),...
                'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.fits']);
        end
    end
    
%     norm_conv(t) = norm(xsol(:) - xstar(:));
%     fprintf('convergence metric: %f\n', norm_conv(t));
    
    if (param.verbose <= 0.5)
        fprintf('.\n');fprintf('\b');
        if mod(t, 50) == 0
            fprintf('\n');
        end
    end
%     if (param.verbose >= 0.5)
%         L1_v(t) = sum(cell2mat(norm1));
%         L2_v(t) = norm(cell2mat(norm2));
%         for q = 1:R
%             L2_vp(t, q) = norm2{q};
%         end
        
%         try 
%             snr_v(t) = 20*log10(norm(param.im(:))/norm(param.im(:) - xsol(:)));
%         catch
%             snr_v(t) = 0;
%         end
%     end
    
    if (param.step_flag && (norm(cell2mat(norm2)) <= param.adapt_eps_tol_out*norm(cell2mat(epsilon))) && ...
       rel_fval(t) < param.reweight_rel_obj && t > 300)
        reweight_steps = [t: param.reweight_step_size :param.max_iter+(2*param.reweight_step_size)];
        param.step_flag = 0;
    end
    
    
%     if (param.use_reweight_steps && t == reweight_steps(reweight_step_count+1) && t < param.reweight_max_reweight_itr) || ...
%         (param.use_reweight_eps && rel_sol_norm_change < param.reweight_rel_obj && ...
%         residual_check_c <= param.adapt_eps_tol_out*epsilon_check_c && residual_check_a <= param.adapt_eps_tol_out*epsilon_check_a && ...
%         param.reweight_min_steps_rel_obj < t - reweight_last_step_iter && t < param.reweight_max_reweight_itr) || ...
%         (param.use_reweight_eps && t - reweight_last_step_iter > param.rw_tol)
    
    if (param.use_reweight_steps && t == reweight_steps(reweight_step_count+1) && t < param.reweight_max_reweight_itr) || ...
        (param.use_reweight_eps && rel_sol_norm_change < param.reweight_rel_obj && ...
        residual_check <= param.adapt_eps_tol_out*epsilon_check && ...
        param.reweight_min_steps_rel_obj < t - reweight_last_step_iter && t < param.reweight_max_reweight_itr) || ...
        (param.use_reweight_eps && t - reweight_last_step_iter > param.rw_tol)
    
        fprintf('\n\n\n\n\n\n\n Performed reweight no %d \n\n\n\n\n', reweight_step_count);
        if (param.verbose >= 1)
            fprintf('Iter %i\n',t);
            fprintf(' L1 norm                       = %e\n', sum(cell2mat(norm1)));
            fprintf(' Residual                      = %e\n', norm(cell2mat(norm2)));
            fprintf(' Epsilon                       = %e\n', norm(cell2mat(epsilon)));
            fprintf(' Relative solution norm change = %e\n\n', rel_sol_norm_change);

            if (param.verbose >= 2)
                for q = 1:R
                    fprintf('   Residual %i                     = %e\n', q, norm2{q});
                    fprintf('   Epsilon %i             = %e\n', q, epsilon{q});
                end
            end
        end
    
%         weights_mat = [weights{:}];
%         weights_mat = weights_mat(:);
%         sigma1 = 1/op_norm(@(x) weights_mat .* Psitw(x), @(x) Psiw(weights_mat .* x), [Ny, Nx], 1e-8, 200, 0);
        
        if reweight_step_count == 0 || reweight_step_count == 1 || ~mod(reweight_step_count,reweight_monitor_step)
%             % Calculate residual image
%             g_f2 = zeros(No, 1);
%             for q = 1 : R
%                 if reduction_version == 1   % H is self-adjoint in this case
%                     tmp = zeros(size(Wm{q}));
%                     tmp(Wm{q}) = T{q} .* res{q};
%                     if flagW
%                         g_f2(W{q}) = g_f2(W{q}) + H{q} * A(real(IFT2(reshape(tmp, Ny, Nx))));
%                     else
%                         g_f2 = g_f2 + H{q} * A(real(IFT2(reshape(tmp, Ny, Nx))));
%                     end
%                 elseif reduction_version == 2
%                     if flagW
%                         g_f2(W{q}) = g_f2(W{q}) + H{q}' * (T{q} .* res{q});
%                     else
%                         g_f2 = g_f2 + H{q}' * (T{q} .* res{q});
%                     end
%                 end
%             end
%             res_im = real(At(g_f2));
            
            % parallel for all bases
            parfor k = 1:P
                d_val = abs(Psit{k}(xsol));
%                 res_val = Psit{k}(res_im/peak_psf);
%                 noise_val = std(res_val);
%                 weights{k} = noise_val * max(reweight_alpha, 1) ./ (noise_val * max(reweight_alpha, 1) + d_val);
                weights{k} = reweight_alpha ./ (reweight_alpha + d_val);
%                 weights{k} = max(reweight_alpha, noise_val) ./ (max(reweight_alpha, noise_val) + d_val);
%                 kstd = 3*1.4826*mad(res_val);
%                 weights{k} = 1 ./ (1e-10 + d_val/(reweight_alpha * gamma));
%                 weights{k}(d_val > max(d_val) * reweight_abs_of_max) = 0;
            end

            reweight_alpha = reweight_alpha_ff * reweight_alpha;
            param.reweight_alpha = reweight_alpha;
            param.init_reweight_step_count = reweight_step_count+1;
            param.init_reweight_last_iter_step = t;
            param.init_t_start = t+1;
            
            if wterm
                fitswrite(xsol, ['results/xsol_ch', num2str(param.chInd),'_prec_it', num2str(t), '_reweight', num2str(reweight_step_count), '_rw_alpha', num2str(reweight_alpha), '_gamma', num2str(gamma)...
                    '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '.fits']);
                fitswrite(res_im, ['results/res_ch', num2str(param.chInd),'_prec_it', num2str(t), '_reweight', num2str(reweight_step_count), '_rw_alpha', num2str(reweight_alpha), '_gamma', num2str(gamma)...
                    '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '.fits']);
            else
                fitswrite(xsol, ['results/xsol_ch', num2str(param.chInd),'_prec_it', num2str(t), '_reweight', num2str(reweight_step_count), '_rw_alpha', num2str(reweight_alpha), '_gamma', num2str(gamma)...
                    '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.fits']);
                fitswrite(res_im, ['results/res_ch', num2str(param.chInd),'_prec_it', num2str(t), '_reweight', num2str(reweight_step_count), '_rw_alpha', num2str(reweight_alpha), '_gamma', num2str(gamma)...
                    '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.fits']);
            end
           
            if wterm
                m = matfile(['./results/SARA_ch', num2str(param.chInd),'_prec_reweight', num2str(reweight_step_count), '_rw_alpha', num2str(reweight_alpha), '_gamma', num2str(gamma)...
                    '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '.mat'], ...
                    'Writable', true);
            else
                m = matfile(['./results/SARA_ch', num2str(param.chInd),'_prec_reweight', num2str(reweight_step_count), '_rw_alpha', num2str(reweight_alpha), '_gamma', num2str(gamma)...
                    '_', num2str(realdatablocks), 'b_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '.mat'], ...
                    'Writable', true);
            end
            m.param = param;
            m.xsol = xsol;
            m.epsilon = epsilon;
            m.weights = weights;
            m.norm1 = norm1;
            m.u1 = u1;
            m.v1 = v1;
            m.norm2 = norm2;
            m.r2 = r2;
            m.u2 = u2;
            m.v2 = v2;
            m.proj = proj;
            m.t_block = t_block;
            m.count_eps_update_down = count_eps_update_down;
            m.count_eps_update_up = count_eps_update_up;
            m.rel_fval = rel_fval;
            m.end_iter = end_iter;
            m.t_l11 = t_l11;
            m.t_master = t_master;
            m.t_data = t_data;
            m.g1 = g1;
            m.g2 = g2;
            m.flag = flag;  
            m.reweight_alpha = reweight_alpha;
%             m.peak_psf = peak_psf;
        end
       
        reweight_last_step_iter = t;
        reweight_step_count = reweight_step_count + 1;
        
        if (reweight_step_count > param.total_reweights)
            param.reweight_max_reweight_itr = t+1;
            fprintf('\n\n No more reweights \n\n');
            break;
        end
    end

    %% global stopping criteria
%     if rel_sol_norm_change < param.rel_obj && ...
%             ((param.global_stop_bound && norm(cell2mat(norm2)) <= norm(cell2mat(epsilon))) || ...
%             (~param.global_stop_bound && prod(cell2mat(norm2) <= cell2mat(epsilon))))
%         flag = 1;
%         break;
%     end
    
%     if (t > 1 && rel_sol_norm_change < param.rel_obj && reweight_step_count > param.total_reweights && ...
%             residual_check_c <= param.adapt_eps_tol_out*epsilon_check_c && residual_check_a <= param.adapt_eps_tol_out*epsilon_check_a)
%         flag = 1;
%         break;
%     end
    if (t > 1 && rel_sol_norm_change < param.rel_obj && reweight_step_count > param.total_reweights && ...
            residual_check <= param.adapt_eps_tol_out*epsilon_check)
        flag = 1;
        break;
     end
end
toc(start_loop)

% final log
if (param.verbose > 0)
    if (flag == 1)
        fprintf('\nSolution found\n');
        fprintf(' L1 norm                       = %e\n', sum(cell2mat(norm1)));
        fprintf(' Residual                      = %e\n', norm(cell2mat(norm2)));
        fprintf(' Epsilon                       = %e\n', norm(cell2mat(epsilon)));
        fprintf(' Relative solution norm change = %e\n\n', rel_sol_norm_change);
        
        for q = 1:R
            fprintf('   Residual %i                     = %e\n', q, norm2{q});
            fprintf('   Epsilon %i             = %e\n', q, epsilon{q});
        end
    else
        fprintf('\nMaximum number of iterations reached\n');
        fprintf(' L1 norm                       = %e\n', sum(cell2mat(norm1)));
        fprintf(' Residual                      = %e\n', norm(cell2mat(norm2)));
        fprintf(' Epsilon                       = %e\n', norm(cell2mat(epsilon)));
        fprintf(' Relative solution norm change = %e\n\n', rel_sol_norm_change);
        
        for q = 1:R
            fprintf('   Residual %i                     = %e\n', q, norm2{q});
            fprintf('   Epsilon %i             = %e\n', q, epsilon{q});
        end
    end
end
fprintf('\n');

xsol = hardt(xsol);

% trim the log vectors to the size of the actual iterations performed
if (param.verbose >= 0.5)
%     L1_v = L1_v(1:t);
%     L2_v = L2_v(1:t);
%     L2_vp = L2_vp(1:t, :);
    rel_fval = rel_fval(1:t);
%     snr_v = snr_v(1:t);
end

end

function [v1_, u1_, l1_, t_l11_] = run_par_waverec(v1_, Psit, Psi, prev_xsol, gamma, weights_, sigma1, lambda1)
    tw = tic;
    r1_ = v1_ + Psit(prev_xsol);
    vy1_ =  r1_ - sign(r1_) .* max(abs(r1_) - gamma * weights_/ sigma1, 0);
    v1_ = v1_ + lambda1 * (vy1_ - v1_);
    u1_ = Psi(v1_);
    t_l11_ = toc(tw);
    % local L11 norm of current solution
    l1_ = norm(r1_(:),1);
end

