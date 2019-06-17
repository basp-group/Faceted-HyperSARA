 function p = prox_negative_root(x, gamma, q)
%function p = prox_negative_root(x, gamma, q)
%
% This procedure computes the proximity operator of the function:
%
%                  f(x) = -gamma * x^(1/q)     with q >= 1
%
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  q     - not less than 1, scalar
% 
%  DEPENDENCIES
% ==============
%  solver3.m - located in the folder 'utils'
%  newton.m  - located in the folder 'utils'
%
%  REMARKS
% =========
% For q ~= {1,2}, the computation is based on Newton's method. Depending 
% on the input values, you may experience slowness or numerical issues.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (27-04-2017)
% Author  : Emilie Chouzenoux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2017
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check input
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x''')
end
if any( q(:) < 1 ) || ~isscalar(q)
    error('''q'' must be a real number not less than 1')
end
%-----%


% compute the prox
if q == 1
    
    p = max(0, x+gamma);
    
elseif q == 2
    
    p = solver3(1, -x, 0, -gamma/2);
    p = sqrt(p);
    
else
     
    % prepare the Newton's method
    fun = @(t) t.^(2*q-1) - x .* t.^(q-1) - gamma/q;
    der = @(t) (2*q-1) * t.^(2*q-2) - (q-1) * x .* t.^(q-2);
    
    % initialize the solution
    p_init = max(1, x);
    
    % use the Newton's method
    p = newton(fun, der, p_init);
    
    % compute the n-th root
    p = nthroot(p, q);

end