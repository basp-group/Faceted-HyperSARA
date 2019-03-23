 function p = prox_max(x, gamma, dir)
%function p = prox_max(x, gamma, dir)
%
% This procedure computes the proximity operator of the function:
%
%                   f(x) = gamma * max(x)
%
% When the input 'x' is an array, the computation can vary as follows:
%  - dir = 0 --> 'x' is processed as a single vector [DEFAULT]
%  - dir > 0 --> 'x' is processed block-wise along the specified direction
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array compatible with the blocks of 'x'
%  dir   - integer, direction of block-wise processing

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
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && (isempty(dir) || any(size(gamma)~=sz))
    error('''gamma'' must be positive and either scalar or compatible with the blocks of ''x''')
end
%------%


%
% 0. Linearize
%
sz = size(x);
if isempty(dir)
    x     = x(:);
    dir = 1;
end


%
% 1. Order the column elements
%
s = sort(x, dir, 'descend');
 

%
% 2. Compute the partial sums: c(j) = ( sum(s(1:j)) - gamma ) / j
%
c = ( cumsum(s,dir) - gamma ) ./ cumsum( ones(size(s)), dir );


%
% 3. Find the index: n = max{ j \in {1,...,B} : s(j) > c(j) }
%
mask = s > c;
mask = mask .* reshape( 1:numel(mask), size(mask) );
n = max(mask, [], dir);


%
% 4. Compute the prox
%
p = bsxfun( @min, x, c(n) );


%
% 5. revert back
%
p = reshape(p, sz);