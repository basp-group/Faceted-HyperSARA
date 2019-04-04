 function p = project_hyperplane(x, eta, w, dir)
%function p = project_hyperplane(x, eta, w, dir)
%
% This procedure computes the projection onto the constraint set:
%
%                           w'x = eta
%
% When the input 'x' is an array, the computation can vary as follows:
%  - dir = 0 --> 'x' is processed as a single vector [DEFAULT]
%  - dir > 0 --> 'x' is processed block-wise along the specified direction
%
%  INPUTS
% ========
%  x   - ND array
%  eta - scalar or ND array compatible with the blocks of 'x'
%  w   - ND array of the same size as 'x' [OPTIONAL]
%  dir - integer, direction of block-wise processing [OPTIONAL]

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
if nargin < 3 || isempty(w)
    w = ones(size(x));
end
if nargin < 4 || (~isempty(dir) && dir == 0)
    dir = [];
end

% check input
sz = size(x); sz(dir) = 1;
if ~isscalar(eta) && any(size(eta) ~= sz)
    error('''eta'' must be either scalar or compatible with the blocks of ''x''')
end
if any( size(w) ~= size(x) )
    error('''w'' must be the same size as ''x''')
end
%-----%


% linearize
sz = size(x);
if isempty(dir)
    x   = x(:);
    w   = w(:);
end

% compute the projection
mu = ( eta - sum(w.*x,dir) ) ./ sum(w.^2,dir);
p = x + bsxfun(@times, mu, w);

% revert back
p = reshape(p, sz);