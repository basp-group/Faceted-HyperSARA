 function [p,q] = prox_diff(x, y, gamma)
%function [p,q] = prox_diff(x, y, gamma)
%
% This procedure computes the proximity operator of the function:
%
%                 / gamma * |x-y|    if x >= 0 and y >= 0
%        D(x,y) = |                                 
%                 \ + inf            otherwise
%
% When the inputs are arrays, the outputs are computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array with the same size as 'y'
%  y     - ND array with the same size as 'x'
%  gamma - positive, scalar or ND array

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
    error('''gamma'' must be positive and either scalar or the same size as ''x''')
end
%-----%


% 4th branch
p = zeros( size(x) );
q = p;

% 3rd branch
mask = (y > gamma) & (x <= -gamma);
yy = y - gamma;
q(mask) = yy(mask);

% 2nd branch
mask = (x > gamma) & (y <= -gamma);
xx = x - gamma;
p(mask) = xx(mask);

% 1st branch
t = sign(x-y) .* max(0, abs(x-y) - 2*gamma);
mask = abs(t) < x+y;
xy = x(mask) + y(mask);
tt = t(mask);
p(mask) = 0.5 * (xy + tt);
q(mask) = 0.5 * (xy - tt);