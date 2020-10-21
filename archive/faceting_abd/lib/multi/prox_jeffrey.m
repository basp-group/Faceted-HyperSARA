 function [p,q] = prox_jeffrey(x, y, gamma)
%function [p,q] = prox_jeffrey(x, y, gamma)
%
% This procedure computes the proximity operator of the function:
%
%             / gamma * (x-y) * log(x/y)  if x > 0 and y > 0
%    D(x,y) = | 0                         if x = y = 0
%             \ +inf                      otherwise
%
% When the inputs are arrays, the outputs are computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array with the same size as 'y'
%  y     - ND array with the same size as 'x'
%  gamma - positive, scalar or ND array
% 
%  DEPENDENCIES
% ==============
%  newton.m  - located in the folder 'utils'
%  Lambert_W.m  - located in the folder 'utils'

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
if any( size(x) ~= size(y) )
    error('''x'' and ''y'' must be the same size')
end
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x'' and ''y''')
end
%-----%


% 2nd branch
p = zeros(size(x));
q = p;

% branch selection
chi_minus = Wexp(1 - x./gamma);
chi_plus  = Wexp(1 - y./gamma);
mask = chi_minus .* chi_plus < 1;
xx = x(mask);
yy = y(mask);
if isscalar(gamma)
    gg = gamma;
else
    gg = gamma(mask);
end

% newton's method
fun = @(t) (t+1) .* log(t) - 1./t    + t.^2 + (xx./gg - 1) .* t + 1 - yy./gg;
der = @(t)  log(t) + 1./t  + 1./t.^2 + 2*t  +  xx./gg;

% root finding
t = chi_minus(mask) + 1e-7;
t = newton(fun, der, t);

% 1st branch
p(mask) = xx + gg .* (log(t) +    t - 1);
q(mask) = yy - gg .* (log(t) - 1./t + 1);



function z = Wexp(u)

z = Lambert_W(exp(u));

% check for big values
mask = exp(u) > realmax / 1e4;

% approximation
v = u(mask);
z(mask) = v - log(v)  .* (1 - 0.5./v);