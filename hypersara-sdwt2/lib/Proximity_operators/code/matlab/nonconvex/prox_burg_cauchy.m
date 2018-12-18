 function p = prox_burg_cauchy(x, gamma, w, d)
%function p = prox_burg_cauchy(x, gamma, w, d)
%
% The function computes the proximity operator of the function:
%
%              f(x) = gamma * ( -log(x) + w*log(d+x^2) )
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
if any( w(:) <= 0 ) || ~isscalar(w) && any(size(w) ~= size(x))
    error('''w'' must be positive and either scalar or the same size as ''x''')
end
if any( d(:) <= 0 ) || ~isscalar(d) && any(size(d) ~= size(x))
    error('''d'' must be positive and either scalar or the same size as ''x''')
end
%-----%


% prepare the Newton's method
fun = @(t)   t.^4 -   x .* t.^3 +   (d+2*gamma.*w-gamma) .* t.^2 - d.*x .* t - d.*gamma;
der = @(t) 4*t.^3 - 3*x .* t.^2 + 2*(d+2*gamma.*w-gamma) .* t    - d.*x;
    
% initialize the solution
p_init = max(1, x);
    
% use the Newton's method
p = newton(fun, der, p_init);