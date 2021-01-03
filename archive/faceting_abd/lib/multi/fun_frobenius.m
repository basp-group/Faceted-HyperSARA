 function p = fun_frobenius(x, gamma)
%function p = fun_frobenius(x, gamma)
%
% This procedure evaluates the function:
%
%                    f(x) = gamma * ||x||_F
% 
% It is assumed that the matrices are stored along the dimensions 1 and 2.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array compatible with the size of 'x'

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
x = reshape(x, [sz(1)*sz(2) 1 sz(3:end)]);

% evaluate the function
xx = sqrt( sum(x.^2) );
p = sum( gamma(:) .* xx(:) );