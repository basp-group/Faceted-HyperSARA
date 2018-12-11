 function y = prox_frobenius(x, gamma)
%function y = prox_frobenius(x, gamma)
%
% This procedure computes the proximity operator of the function:
%
%                    f(x) = gamma * ||x||_F
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
%  prox_L2.m  - located in the folder 'multi'

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


% reshape the matrices as columns
sz = size(x);
s = reshape(x, [sz(1)*sz(2) 1 sz(3:end)]);

% compute the vectorial prox
g = prox_L2(s, gamma, 1);

% reshape the columns as matrices
y = reshape(g, sz);