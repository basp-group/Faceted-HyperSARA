 function p = prox_moreau(x, gamma, b, prox_phi, varargin)
%function p = prox_moreau(x, gamma, b, prox_phi, varargin)
%
% This procedure computes the proximity operator of a function:
% 
%          f(x) = gamma * inf_y  phi(y) + (y-x)^2 / (2*b)
%
% When the input 'x' is an array, the expression '(y-x)^2 / (2*b)' is meant
% element-wise. However, from a theoretical standpoint, the input 'b' must
% be scalar if 'prox_phi' represents a non-separable proximity operator.
%
%  INPUTS
% ========
% x        - ND array
% gamma    - positive, scalar or ND array with the same size as 'x'
% b        - positive, scalar or ND array with the same size as 'x'
% prox_phi - function handle with two arguments at least
% varargin - additional parameters for the function 'prox_phi' [OPTIONAL]

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



% check inputs
if any( b(:) < 0 ) || ~isscalar(b) && any(size(b) ~= size(x))
    error('''b'' must be scalar or the same size as ''x''');
end
%-----%


% compute the prox
p = prox_phi(x, b+gamma, varargin{:});
p = (b .* x + gamma .* p) ./ (b+gamma);