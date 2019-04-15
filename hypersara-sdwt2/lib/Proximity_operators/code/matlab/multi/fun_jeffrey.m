 function p = fun_jeffrey(x, y, gamma)
%function p = fun_jeffrey(x, y, gamma)
%
% This procedure evaluates the function:
%
%             / gamma * (x-y) * log(x/y)  if x > 0 and y > 0
%    D(x,y) = | 0                         if x = y = 0
%             \ +inf                      otherwise
%
% When the inputs are arrays, the output is the element-wise sum.
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
    error('''gamma'' must be positive and either scalar or the same size as ''x'' and ''y''')
end
%-----%


if any( x(:)<0 | y(:)<0 )
    p = inf;
else
    p = gamma .* (x-y) .* log(x./y);
    p(x == 0 & y == 0) = 0;
    p = sum(p(:));
end