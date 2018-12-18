 function p = prox_nuclear(x, gamma)
%function p = prox_nuclear(x, gamma)
%
% This procedure computes the proximity operator of the function:
%
%                    f(x) = gamma * ||x||_N
% 
% It is assumed that the matrices are stored along the dimensions 1 and 2.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array compatible with the size of 'x'
% 
%  DEPENDENCIES
% ==============
%  prox_abs.m     - located in the folder 'scalar'
%  prox_svd.m     - located in the folder 'utils'
%  sv_dec.cpp     - [MEX-FILE] located in the folder 'utils'
%  sv_dec_fat.cpp - [MEX-FILE] located in the folder 'utils'
%  sv_rec.cpp     - [MEX-FILE] located in the folder 'utils'

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


p = prox_svd(x, gamma, @prox_abs);