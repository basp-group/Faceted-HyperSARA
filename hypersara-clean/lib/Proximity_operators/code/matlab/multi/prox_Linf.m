 function p = prox_Linf(x, gamma, dir)
%function p = prox_Linf(x, gamma, dir)
%
% This procedure computes the proximity operator of the function:
%
%                   f(x) = gamma * ||x||_inf
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
% 
%  DEPENDENCIES
% ==============
%  prox_max.m  - located in the folder 'multi'

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
if nargin < 3, dir = []; end
%-----%

p = sign(x) .* prox_max( abs(x), gamma, dir );