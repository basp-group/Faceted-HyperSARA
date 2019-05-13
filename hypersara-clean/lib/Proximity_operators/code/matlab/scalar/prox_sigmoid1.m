 function p = prox_sigmoid1(x )
%function p = prox_sigmoid1(x )
%
% This procedure computes the proximity operator of the function:
%
%                       f(x) = -\sqrt{1 - x^2} - 1/2 x^2 
%
% When the input 'x' is an array, the output 'p' is computed ...
%
%  INPUTS
% ========
%  x     - ND array 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (dd-mm-yyyy)
% Author  : xxx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2017
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 


% compute the prox
p = x./sqrt(1 +x.^2);