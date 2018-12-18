 function p = project_simplex_ext(x, eta, dir)
%function p = project_simplex_ext(x, eta, dir)
%
% This procedure computes the projection onto the constraint set:
%
%                  x => 0   AND   1'x <= eta
%
% When the input 'x' is an array, the computation can vary as follows:
%  - dir = 0 --> 'x' is processed as a single vector [DEFAULT]
%  - dir > 0 --> 'x' is processed block-wise along the specified direction
%
%  INPUTS
% ========
%  x   - ND array
%  eta - positive, scalar or ND array compatible with the blocks of 'x'
%  dir - integer, direction of block-wise processing
% 
%  DEPENDENCIES
% ==============
%  project_simplex.m  - located in the folder 'indicator'
%  prox_max.m         - located in the folder 'multi'

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
if nargin < 3 || (~isempty(dir) && dir == 0)
    dir = [];
end

% check input
sz = size(x); sz(dir) = 1;
if any( eta(:) <= 0 ) || ~isscalar(eta) && any(size(eta) ~= sz)
    error('''eta'' must be positive and either scalar or compatible with the blocks of ''x''')
end
%-----%

% reshape
sz = size(x);
if isempty(dir)
    x = x(:);
end

% pre-project onto the positive orthant
p = max(x,0);

% find the blocks still outside the simplex
mask = sum(p,dir) > eta;

% project onto the simplex
pp = project_simplex(p, eta, dir);

% copy the projected blocks
rep = 1 + (size(x)-1) .* (1:ndims(x) == dir);
mask = repmat(mask, rep);
p(mask) = pp(mask);

% revert back
if isempty(dir)
    p = reshape(p, sz);
end