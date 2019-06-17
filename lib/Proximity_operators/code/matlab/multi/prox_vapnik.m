 function p = fun_vapnik(x, gamma, w, dir)
%function p = fun_vapnik(x, gamma, w, dir)
%
% This procedure computes the proximity operator of the function:
%
%                f(x) = gamma * max(||x||_2 - w)             with w > 0
%
% When the input 'x' is an array, the computation can vary as follows:
%  - dir = 0 --> 'x' is processed as a single vector [DEFAULT]
%  - dir > 0 --> 'x' is processed block-wise along the specified direction
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array compatible with the blocks of 'x'
%  w     - positive, scalar or ND array with the same size as 'x'
%  dir   - integer, direction of block-wise processing

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
sz = size(x); sz(dir) = 1;
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && (isempty(dir) || any(size(gamma)~=sz))
    error('''gamma'' must be positive and either scalar or compatible with the blocks of ''x''')
end
if any( w(:) <= 0 ) || ~isscalar(w) && (isempty(dir) || any(size(w)~=sz))
    error('''w'' must be positive and either scalar or compatible with the blocks of ''x''')
end
%------%


% linearize
sz = size(x);
if isempty(dir)
    x = x(:);
end
    
% preliminaries
norm_x = sqrt( sum(x.^2, dir) );

% 3rd branch
p = bsxfun(@times, x, 1-gamma./norm_x);

% 2nd branch
mask = (w < norm_x) & (norm_x <= w+gamma);
p2 = bsxfun(@times, x, w ./ norm_x);
p(mask) = p2(mask);

% 1st branch
mask = norm_x <= w;
p(mask) = x(mask);

% revert back
p = reshape(p, sz);