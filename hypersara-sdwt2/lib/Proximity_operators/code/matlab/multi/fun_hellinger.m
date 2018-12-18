 function p = fun_hellinger(x, y, gamma)
%function p = fun_hellinger(x, y, gamma)
%
% This procedure evaluates the function:
%
%            / gamma * ( sqrt(x) - sqrt(y) )^2   if x >= 0 and y >= 0
%   D(x,y) = |
%            \ +inf                              otherwise
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



if any( x(:)<0 | y(:)<0 )
    p = inf;
else
    p = gamma .* ( sqrt(x) - sqrt(y) ).^2;
    p = sum(p(:));
end