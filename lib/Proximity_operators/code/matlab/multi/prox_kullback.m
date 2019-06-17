 function [p,q] = prox_kullback(x, y, gamma, kappa)
%function [p,q] = prox_kullback(x, y, gamma, kappa)
%
% This procedure computes the proximity operator of the function:
%
%           / gamma * ( x * log(x/y) + k * (y - x) )   if x > 0 and y > 0
%  D(x,y) = | 0                                        if x = 0 and y >=0
%           \ +inf                                     otherwise
%
% When the inputs are arrays, the outputs are computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array with the same size as 'y'
%  y     - ND array with the same size as 'x'
%  gamma - positive, scalar or ND array
%  kappa - any real, scalar or ND array [DEFAULT: kappa=1]
% 
%  DEPENDENCIES
% ==============
%  newton.m  - located in the folder 'utils'

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


% default input
if nargin < 4 || isempty(kappa)
    kappa = 1;
end

% check input
if any( size(x) ~= size(y) )
    error('''x'' and ''y'' must be the same size')
end
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x'' and ''y''')
end
if ~isscalar(kappa) && any(size(kappa) ~= size(x))
    error('''kappa'' must be either scalar or the same size as ''x'' and ''y''')
end
%-----%


% input correction
x = x + gamma * kappa - gamma;
y = y - gamma * kappa + gamma;

% 2nd branch
p = zeros(size(x));
q = p;

% branch selection
mask1 = (y >= gamma) | ( x./gamma > log(1 - y./gamma) );
mask2 = exp(-x./gamma) > 1/eps;

% approximation when exp(-x/gamma) == inf
p(mask1 & mask2) = 0;
yy = y - gamma;
q(mask1 & mask2) = yy(mask1 & mask2);

% branch selection
mask = mask1 & ~mask2;
xx = x(mask);
yy = y(mask);
if isscalar(gamma)
    gg = gamma;
else
    gg = gamma(mask);
end

% newton's method
fun = @(t) t .* log(t) + xx .* t ./ gg - 1./t    + 1 - yy./gg;
der = @(t) 1 +  log(t) + xx      ./ gg + 1./t.^2;

% root finding
t = exp(-xx./gg) + eps;
t = newton(fun, der, t);

% 1st branch
p(mask) = xx + gg .* log(t);
q(mask) = yy + gg .* (1./t - 1);