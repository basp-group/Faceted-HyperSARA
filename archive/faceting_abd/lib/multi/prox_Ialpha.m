 function [p,q] = prox_Ialpha(x, y, gamma, kappa)
%function [p,q] = prox_Ialpha(x, y, gamma, kappa)
%
% This procedure computes the proximity operator of the function:
%
%          / gamma * ( -sqrt(x * y) + k/2 * (x + y) )  if x >= 0 and y >= 0
% D(x,y) = |
%          \ +inf
%
% When the inputs are arrays, the outputs are computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array with the same size as 'y'
%  y     - ND array with the same size as 'x'
%  gamma - positive, scalar or ND array
%  kappa - any real, scalar or ND array [DEFAULT: kappa=0]

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
    kappa = 0;
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
x = x + (1-kappa) .* gamma/2;
y = y + (1-kappa) .* gamma/2;

% 2nd branch
p = zeros( size(x) );
q = p;

% branch selection
mask = (x >= gamma/2) | ( (1 - 2*y ./ gamma) < 1 ./ (1 - 2*x./gamma) );
xx = x(mask);
yy = y(mask);
if isscalar(gamma)
    gg = gamma;
else
    gg = gamma(mask);
end

% newton's method
a = xx./gg - 0.5;
b = 0.5 - yy./gg;
fun = @(t) 0.5 * t.^4 +   a .* t.^3 + b .* t - 0.5;
der = @(t)   2 * t.^3 + 3*a .* t.^2 + b;

% root finding
rho = 100 + max(0, 1-2*xx./gg);
rho = newton(fun, der, rho);

% special cases
rho(a==0 & b==0) = 1;

% 1st branch
p(mask) = xx + gg/2 .* (   rho - 1);
q(mask) = yy + gg/2 .* (1./rho - 1);