 function p = prox_entropy_logsum(x, gamma, w, d)
%function p = prox_entropy_logsum(x, gamma, w, d)
%
% The function computes the proximity operator of the function:
%
%              / gamma * ( x * log(x) + w*log(d+x) )  if x > 0
%       f(x) = |
%              \ +inf
%
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  w     - positive, scalar or ND array with the same size as 'x'
%  d     - positive, scalar or ND array with the same size as 'x'
% 
%  DEPENDENCIES
% ==============
%  newton.m - located in the folder 'utils'
%
%  REMARKS
% =========
% For some inputs, we need to find multiple zeros of a non-linear equation.
% In these cases, the Newton's method just picks one of them. Hence, there 
% is no guarantee that we pick the zero that leads to the correct output !

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
if any( w(:) <= 0 ) || ~isscalar(w) && any(size(w) ~= size(x))
    error('''w'' must be positive and either scalar or the same size as ''x''')
end
if any( d(:) <= 0 ) || ~isscalar(d) && any(size(d) ~= size(x))
    error('''d'' must be positive and either scalar or the same size as ''x''')
end
%-----%


% prepare the Newton's method
fun = @(t) t.^2 + (d-x+gamma).*t + gamma.*d.*log(t) + gamma.*t.*log(t) + d.*(gamma-x) + w.*gamma;
der = @(t)  2*t + (d-x+gamma)    + gamma.*d./t      + gamma.*(log(t)+1);
    
% initialize the solution
p_init = max(1, x);
    
% use the Newton's method
p = newton(fun, der, p_init);