 function [p,q] = prox_perspective_square(y, xi, gamma, dir)
%function [p,q] = prox_perspective_square(y, xi, gamma, dir)
%
% This procedure computes the proximity operator of the function
%
%                 / gamma * ||y||_2^2 / xi   if xi > 0
%       f(y,xi) = | 0                        if y = 0 and xi = 0
%                 \ +inf                     otherwise
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

% 3rd branch
yy = sqrt( sum(y.^2,dir) );
t = solver3(gamma, 0, 4*xi+8*gamma, - 8*yy);
pp = gamma .* t   ./ yy;
qq = gamma .* t.^2 / 4;

% 1st branch
mask = yy.^2 <= -4*gamma.*xi;
pp(mask) = 1;
qq(mask) = -xi(mask);

% 2nd branch
mask = yy==0 & xi > 0;
pp(mask) = 1;
qq(mask) = 0;

% compute the prox
p =  y - bsxfun(@times, y, pp);
q = xi + qq;

% revert back
p = reshape(p, sz);