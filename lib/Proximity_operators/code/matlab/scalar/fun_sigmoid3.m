 function p = fun_sigmoid3(x)
%function p = fun_sigmoid3(x)
%
% This procedure evaluates the function
%
% f(x) =  [-2/pi log(cos(pi/2 x)) -1/2 x^2]
%
% When the input 'x' is an array, the output 'p' is ...
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

% evaluate the function
xx =x(:);
p = sum(-2/pi.*log(cos(pi/2.*xx)) -1/2.*xx.^2);
p(abs(xx) > 1) = Inf;