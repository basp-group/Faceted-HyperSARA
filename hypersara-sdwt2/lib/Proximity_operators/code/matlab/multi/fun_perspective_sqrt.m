 function p = fun_perspective_sqrt(y, xi, gamma, dir)
%function p = fun_perspective_sqrt(y, xi, gamma, dir)
%
% This procedure computes the function
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
if isempty(dir)
    y = y(:);
end

% evaluate the function
yy = sum(y.^2,dir);
if any( xi(:)<0 | yy(:) > 0 )
    p = inf;
else
    p = -gamma .* sqrt(xi.^2 - yy);
    p = sum(p(:));
end