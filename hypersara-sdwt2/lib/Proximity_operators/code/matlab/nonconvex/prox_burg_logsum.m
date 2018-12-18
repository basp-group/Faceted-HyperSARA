 function p = prox_burg_logsum(x, gamma, w, d)
%function p = prox_burg_logsum(x, gamma, w, d)
%
% The function computes the proximity operator of the function:
%
%                  f(x) = gamma * ( -log(x) + w*log(d+x) )
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
%  solver3.m - located in the folder 'utils'

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


% compute the prox
[p,p2,p3] = solver3(1, d-x, w.*gamma-d.*x-gamma, -d.*gamma);

% select the correct solution
fun = @(t) 0.5 * abs(x-t).^2 - gamma .* log(t) + gamma .* w .* log(d+t);
mask = ~imag(p2) & p2 > 0 & fun(p2) < fun(p);
p(mask) = p2(mask);
mask = ~imag(p3) & p3 > 0 & fun(p3) < fun(p);
p(mask) = p3(mask);