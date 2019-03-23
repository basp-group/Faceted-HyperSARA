 function p = fun_vapnik(x, gamma, w, dir)
%function p = fun_vapnik(x, gamma, w, dir)
%
% This procedure evaluates the function:
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
    error('''w'' must be positive and either scalar or the same size as ''x''')
end
%------%


% linearize
if isempty(dir)
    x = x(:);
end
    
% evaluate the function
xx = sqrt( sum(x.^2, dir) );
xx = gamma .* max(xx-w, 0);
p  = sum( xx(:) );