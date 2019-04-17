 function p = prox_cauchy(x, gamma, w)
%function p = prox_cauchy(x, gamma, w)
%
% This procedure computes the proximity operator of the function:
%
%                f(x) = gamma * log(w + x^2)
%
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  w     - positive, scalar or ND array with the same size as 'x'
% 
%  DEPENDENCIES
% ==============
%  solver3.m - located in the folder 'utils'

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
%-----%

 
% solve the cubic equation
[p,p2,p3] = solver3(1, -x, w+2*gamma, -x.*w);

% select the correct solution
fun = @(t) 0.5 * abs(x-t).^2 + gamma .* log(w+t.^2);
mask = ~imag(p2) & fun(p2) < fun(p);
p(mask) = p2(mask);
mask = ~imag(p3) & fun(p3) < fun(p);
p(mask) = p3(mask);