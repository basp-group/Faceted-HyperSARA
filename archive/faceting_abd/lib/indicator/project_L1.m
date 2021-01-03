 function p = project_L1(x, eta, dir)
%function p = project_L1(x, eta, dir)
%
% This procedure computes the projection onto the constraint set:
%
%                         ||x||_1 <= eta
%
% When the input 'x' is an array, the computation can vary as follows:
%  - dir = 0 --> 'x' is processed as a single vector [DEFAULT]
%  - dir > 0 --> 'x' is processed block-wise along the specified direction
%
%  INPUTS
% ========
%  x   - ND array
%  eta - positive, scalar or ND array compatible with the blocks of 'x'
%  dir - integer, direction of block-wise processing [OPTIONAL]
% 
%  DEPENDENCIES
% ==============
%  prox_Linf.m - located in the folder 'multi'
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

p = x - prox_Linf(x, eta, dir);