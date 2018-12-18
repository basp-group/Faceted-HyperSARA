 function p = prox_log_power(x, gamma, w, q)
%function p = prox_log_power(x, gamma, w, q)
%
% This procedure computes the proximity operator of the function:
%
%         / gamma * (-log(x) + w * x^q)  if x > 0
% f(x) = |                                           with w > 0 and q >= 1
%         \ +inf                        otherwise
%
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  w     - positive, scalar or ND array with the same size as 'x'
%  q     - not less than 1, scalar
% 
%  DEPENDENCIES
% ==============
%  solver3.m - located in the folder 'utils'
%  newton.m  - located in the folder 'utils'
%
%  REMARKS
% =========
% Except for q=1,2,3, the computation is based on Newton's method. Hence,
% depending on the inputs, you may experience slowness or numerical issues.

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
if any( w(:) <= 0 ) || ~isscalar(w) && any(size(w) ~= size(x))
    error('''w'' must be positive and either scalar or the same size as ''x''')
end
%-----%


if q == 1
    
    p = ( x - w.*gamma + sqrt((x-w.*gamma).^2 + 4*gamma) ) / 2;
    
elseif q == 2
    
    p = ( x + sqrt(x.^2 + 4*gamma.*(1+2*w.*gamma)) ) / (2 + 4*w.*gamma);
    
elseif q == 3
    
    p = solver3(3*w.*gamma, 1, -x, -gamma);
    
else
    
    % prepare the Newton's method
    fun = @(t) q   * w.*gamma .* t.^q     + t.^2 - x .* t - gamma;
    der = @(t) q^2 * w.*gamma .* t.^(q-1) + t    - x;
    
    % initialize the solution
    p_init = max(1, x);
    
    % use the Newton's method
    p = newton(fun, der, p_init);
    
end