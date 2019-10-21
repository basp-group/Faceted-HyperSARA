function [sol,norm_res] = fb_dr_nnls(y, A, At, H, Sigma, Mask, param)
% fb_nnls - Solve non negative least squares problem.
%
% The problem is solved using the forward-backward algorithm with
% FISTA-like acceleration
% 
% 
% Author: Ming Jiang
%

No = size(Mask, 1);

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end


% Initialization
if isfield(param,'initsol')
    xhat = param.initsol;
    HFx = H * A(xhat);
    res = y - Sigma .* HFx(Mask);
    g2 = zeros(No,1);
    g2(Mask) = g2(Mask) + Sigma .* res;
    grad = At(H * g2);
    prev_obj = 0.5*norm(res(:))^2;
else
    g2 = zeros(No,1);
    g2(Mask) = g2(Mask) + Sigma .* y;
    grad = At(H * g2);
    xhat = zeros(size(grad));
    prev_obj = 0.5*norm(y(:))^2;
end

iter = 1; 
prev_sol = xhat;
told = 1;
beta = 0.5;
qfval = prev_obj;

L2_v = zeros(param.max_iter, 1);

sol_steps = param.sol_steps;
sol_step_count = 1;
sol_v = zeros(length(sol_steps)-1, size(xhat, 1), size(xhat, 2));

snr_v = zeros(param.max_iter, 1);



% Main loop
while 1
    
    %Log
    if param.verbose >= 1
        fprintf('Iteration %i:\n', iter);
    end
    
    %Step size computation    
    HFx = H * A(grad);
    res = Sigma .* HFx(Mask);
    mu = param.beta * norm(grad(:))^2/norm(res(:))^2;
    
    %Gradient descend step
    sol = xhat + mu*grad;
    
    %Projection onto the positive orthant
    sol = real(sol);
    sol(sol<0) = 0;
    
    %Stepsize check
    q = qfval + real((xhat(:)-sol(:))'*grad(:))...
        + 0.5/mu*norm(sol(:)-xhat(:))^2;
    HFx = H * A(sol);
    res = y - Sigma .* HFx(Mask);
    norm_res = norm(res);
    curr_obj = 0.5*norm(res(:))^2;
    
    
    if param.verbose >= 0.5
        L2_v(iter) = norm(res);
    end
    
    while (curr_obj > q)
        %Backtracking rule
        mu = beta*mu;
        %Gradient descend step
        sol = xhat + mu*grad;
    
        %Projection onto the positive orthant
        sol = real(sol);
        sol(sol<0) = 0;
        
        %New stepsize check
        q = qfval + real((sol(:)-xhat(:))'*grad(:))...
            + 0.5/mu*norm(sol(:)-xhat(:))^2;
        HFx = H * A(sol);
        res = y - Sigma .* HFx(Mask);
        norm_res = norm(res);
        curr_obj = 0.5*norm(res(:))^2;  
    end
     
    
    % Global stopping criterion
    rel_obj = abs(curr_obj - prev_obj)/curr_obj;
    if param.verbose >= 1
        fprintf('  Obj = %e, rel_obj = %e\n', ...
            curr_obj, rel_obj);
    end
    if (param.verbose <= 0.5)
        fprintf('.\n');fprintf('\b');
        if mod(iter, 50) == 0
            fprintf('\n');
        end
    end
    if (param.verbose >= 0.5)
        if iter == sol_steps(sol_step_count)
            asol = real(sol);
            asol(asol < 0) = 0;
            sol_v(sol_step_count, :, :) = asol;
            
            sol_step_count = sol_step_count + 1;
        end
        
        asol = real(sol);
        asol(asol < 0) = 0;
        try
            snr_v(iter) = 20*log10(norm(param.im(:))/norm(param.im(:) - asol(:)));
        end
    end
    if (rel_obj < param.rel_obj)
        crit_BPDN = 'REL_OBJ';
        break;
    elseif iter >= param.max_iter
        crit_BPDN = 'MAX_IT';
        break;
    end
    
    % FISTA update
    t = (1+sqrt(1+4*told^2))/2;
    xhat = sol + (told-1)/t * (sol - prev_sol);
    
    % Gradient computation
    HFx = H * A(xhat);
    res = y - Sigma .* HFx(Mask);
    norm_res = norm(res);
    
    g2 = zeros(No,1);
    g2(Mask) = g2(Mask) + Sigma .* res;
    grad = At(H * g2);

    % Update variables
    iter = iter + 1;
    prev_obj = curr_obj;
    prev_sol = sol;
    told = t;
    qfval = 0.5*norm(res(:))^2;
    

end

% Log
if param.verbose >= 0.5
    % L1 norm
    fprintf('\n Solution found:\n');
    fprintf(' Final objective value: %e\n', curr_obj);
    
    % Stopping criterion
    fprintf(' %i iterations\n', iter);
    fprintf(' Stopping criterion: %s \n\n', crit_BPDN);
    
end


% trim the log vectors to the size of the actual iterations performed
if (param.verbose >= 0.5)
    L2_v = L2_v(1:iter-1);
    sol_v = sol_v(1:sol_step_count-1, :, :);
    snr_v = snr_v(1:iter-1);
end

end