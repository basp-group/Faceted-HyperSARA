 function p = prox_inverse(x, gamma, q)
%function p = prox_inverse(x, gamma, q)
%
% This procedure computes the proximity operator of the function:
%
%                  f(x) = gamma / x^q     with q >= 1
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
% For q > 1, the computation is based on Newton's method. Hence, depending 
% on the input values, you may experience slowness or numerical issues.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (27-04-2017)
% Author  : Giovanni Chierchia
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
    
    % find the unique positive root of p^3 - x p^2 - gamma = 0
    p = solver3(1, -x, 0, -gamma);
    
else
     
    % prepare the Newton's method
    fun = @(t) t.^(q+2) - x .* t.^(q+1) - q*gamma;
    der = @(t) (q+2) * t.^(q+1) - (q+1) * x .* t.^q;
    
    % initialize the solution
    p_init = max(1, x);
    
    % use the Newton's method
    p = newton(fun, der, p_init);

end