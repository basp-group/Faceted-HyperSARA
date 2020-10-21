 function [p,q] = prox_renyi(x, y, gamma)
%function [p,q] = prox_renyi(x, y, gamma)
%
% This procedure computes the the proximity operator of the function:
%
%              / gamma * x^2 / y   if x >= 0 and y > 0
%     D(x,y) = | 0                 if x = y = 0
%              \ +inf
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
%  solver3.m  - located in the folder 'utils'

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
p = zeros( size(x) );
q = max(0, y);

% branch selection
mask = (x > 0) & (-4*gamma.*y < x.^2);
xx = x(mask);
yy = y(mask);
if isscalar(gamma)
    gg = gamma;
else
    gg = gamma(mask);
end

% root finding
rho = solver3(1, 0, 2+yy./gg, -xx./gg);

% 1st branch
p(mask) = xx - 2 * gg .* rho;
q(mask) = yy + gg .* rho.^2;