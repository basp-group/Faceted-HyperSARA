 function p = fun_thresh(x, gamma, a, b)
%function p = fun_thresh(x, gamma, a, b)
%
% This procedure evaluates the function:
%
%                  / gamma * a * x   if x < 0
%           f(x) = | 0               if x = 0         with a < b
%                  \ gamma * b * x   otherwise
%
% When the input 'x' is an array, the output 'p' is the element-wise sum.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  a     - scalar or ND array with the same size as 'x'
%  b     - scalar or ND array with the same size as 'x'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (27-04-2017)
% Author  : Emilie Chouzenoux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2017
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% check input
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x''')
end
if ~isscalar(a) && any(size(a) ~= size(x))
    error('''a'' must be either scalar or the same size as ''x''')
end
if ~isscalar(b) && any(size(b) ~= size(x))
    error('''b'' must be either scalar or the same size as ''x''')
end
if ~all( a(:) < b(:) )
    error('''a'' must be less than ''b''')
end
%-----%


% evaluate the function
p = a .* min(0,x) + b .* max(0,x);
p = sum(gamma(:) .* p(:));