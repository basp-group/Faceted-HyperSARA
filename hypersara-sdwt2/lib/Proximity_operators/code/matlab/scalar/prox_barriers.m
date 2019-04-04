 function p = prox_barriers(x, gamma, a, b)
%function p = prox_barriers(x, gamma, a, b)
%
% This procedure computes the proximity operator of the function:
%
%         f(x) = gamma * ( -log(x-a) - log(b-x) )     with a < b
%
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  a     - scalar or ND array with the same size as 'x'
%  b     - scalar or ND array with the same size as 'x'
% 
%  DEPENDENCIES
% ==============
%  solver3.m - located in the folder 'utils'

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


% prepare the cubic equation
c3 = 1;
c2 = -(a+b+x);
c1 = a.*b + (a+b).*x - 2*gamma;
c0 = -a.*b.*x + (a+b).*gamma;

% compute the prox
[~,p,~] = solver3(c3, c2, c1, c0);