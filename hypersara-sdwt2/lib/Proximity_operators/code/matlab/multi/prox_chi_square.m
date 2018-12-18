 function [p,q] = prox_chi_square(x, y, gamma)
%function [p,q] = prox_chi_square(x, y, gamma)
%
% This procedure computes the the proximity operator of the function:
%
%                / gamma * (x - y)^2 / y   if x >= 0 and y > 0
%       D(x,y) = | 0                       if x = y = 0
%                \ +inf                    otherwise
%
% When the inputs are arrays, the outputs are computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array with the same size as 'y'
%  y     - ND array with the same size as 'x'
%  gamma - positive, scalar or ND array
% 
%  DEPENDENCIES
% ==============
%  prox_renyi.m  - located in the folder 'multi'

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


[p,q] = prox_renyi(x+2*gamma, y-gamma, gamma);