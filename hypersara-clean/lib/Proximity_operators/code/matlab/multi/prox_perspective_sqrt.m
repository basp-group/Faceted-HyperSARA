 function [p,q] = prox_perspective_sqrt(y, xi, gamma, dir)
%function [p,q] = prox_perspective_sqrt(y, xi, gamma, dir)
%
% This procedure computes the proximity operator of the function
%
%           / -gamma * sqrt(xi^2 - ||y||_2^2) if xi > 0 and ||y||_2^2 <= xi
% f(y,xi) = | 0                               if y = xi = 0
%           \ +inf                            otherwise
%
% When the input 'y' is an array, the computation can vary as follows:
%
%  - dir = 0 --> 'y' is processed as a single vector [DEFAULT]
%                (in this case, 'xi' must be a scalar) 
%
%  - dir > 0 --> 'y' is processed block-wise along the specified direction
%                (in this case, 'xi' must be singleton along 'dim')
%
%  INPUTS
% ========
%  y     - ND array
%  xi    - ND array compatible with the blocks of 'x'
%  gamma - positive, scalar or ND array with the same size as 'xi'
%  dir   - integer, direction of block-wise processing [DEFAULT: 0]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (26-01-2018)
% Author  : Giovanni Chierchia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2017
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% default inputs
if nargin < 4 || (~isempty(dir) && dir == 0)
    dir = [];
end

% check input
sz = size(y); sz(dir) = 1;
if isempty(dir) && numel(xi) ~= 1 || ~isempty(dir) && any( sz ~= size(xi) )
    error('The input ''xi'' is not compatible with the blocks of ''y''');
end
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && (isempty(dir) || any(size(gamma)~=sz))
    error('''gamma'' must be positive and either scalar or the same size as ''xi''')
end
%-----%


% linearize
sz = size(y);
if isempty(dir)
    y = y(:);
end

% init
p = zeros(size(y));
q = zeros(size(xi));

% condition
yy   = sqrt( sum(y.^2,dir) );
mask = xi + sqrt(gamma.^2 + yy) > 0;
if ~isscalar(gamma)
    gamma = gamma(mask);
end
if ~isscalar(xi)
    xi = xi(mask);
end
yy_mask = yy(mask);

% compute t
fun  = @(t) 2*gamma.*t + xi.*t ./ sqrt(1+t.^2) - yy_mask;
der  = @(t) 2*gamma + (xi.*sqrt(1+t.^2) - xi.*t.^2./sqrt(1+t.^2))./(1+t.^2);
t = newton(fun, der, yy_mask);

% prox 1
coef       = ones(size(yy));
coef(mask) = gamma.*t ./ yy_mask;
p          = y - bsxfun(@times, y, coef);

% prox 2
qq = xi - gamma .* sqrt(1+t.^2);
q(mask) = qq;

% revert back
p = reshape(p, sz);