 function p = fun_burg_logsum(x, gamma, w, d)
%function p = fun_burg_logsum(x, gamma, w, d)
%
% The function evaluates the function:
%
%                  / gamma * ( -log(x) + w*log(d+x) )  if x > 0
%           f(x) = |
%                  \ +inf
%
% When the input 'x' is an array, the output 'p' is the element-wise sum.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  w     - positive, scalar or ND array with the same size as 'x'
%  d     - positive, scalar or ND array with the same size as 'x'

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
if any( w(:) <= 0 ) || ~isscalar(w) && any(size(w) ~= size(x))
    error('''w'' must be positive and either scalar or the same size as ''x''')
end
if any( d(:) <= 0 ) || ~isscalar(d) && any(size(d) ~= size(x))
    error('''d'' must be positive and either scalar or the same size as ''x''')
end
%-----%


% evaluate the function
if any( x(:) <= 0 )
    p = Inf;
else
    p = sum( gamma(:) .* ( -log(x(:)) + w(:) .* log(d(:)+x(:)) ) );
end