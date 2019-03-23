 function p = fun_kullback(x, y, gamma, kappa)
%function p = fun_kullback(x, y, gamma, kappa)
%
% This procedure evaluates the function:
%
%           / gamma * ( x * log(x/y) + k * (y - x) )   if x > 0 and y > 0
%  D(x,y) = | 0                                        if x = 0 and y >=0
%           \ +inf                                     otherwise
%
% When the inputs are arrays, the output is the element-wise sum.
%
%  INPUTS
% ========
%  x     - ND array with the same size as 'y'
%  y     - ND array with the same size as 'x'
%  gamma - positive, scalar or ND array
%  kappa - any real, scalar or ND array [DEFAULT: kappa=1]

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



% default input
if nargin < 4 || isempty(kappa)
    kappa = 1;
end

% check input
if any( size(x) ~= size(y) )
    error('''x'' and ''y'' must be the same size')
end
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x'' and ''y''')
end
if ~isscalar(kappa) && any(size(kappa) ~= size(x))
    error('''kappa'' must be either scalar or the same size as ''x'' and ''y''')
end
%-----%


if any( x(:)<0 | y(:)<0 )
    p = inf;
else
    p = gamma .* ( x .* log(x./y) + k .* (y - x) );
    p(x == 0 & y >= 0) = 0;
    p = sum(p(:));
end
