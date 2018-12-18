 function y = prox_outer(x, gamma, a, b, prox_phi, varargin)
%function y = prox_outer(x, gamma, a, b, prox_phi, varargin)
%
% This procedure computes the proximity operator of a function:
% 
%             f(x) = gamma * ( phi(x) + a * x^2 + b * x )
%
% When the input 'x' is an array, the expression 'a * x^2 + b * x' is meant 
% element-wise. However, from a theoretical standpoint, the input 'a' must 
% be scalar if 'prox_phi' represents a non-separable proximity operator.
%
%  INPUTS
% ========
% x        - ND array
% gamma    - positive, scalar or ND array with the same size as 'x'
% a        - non-negative, scalar or ND array with the same size as 'x' [OPTIONAL]
% b        - scalar or ND array with the same size as 'x' [OPTIONAL]
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


% default inputs
if nargin < 3 || isempty(a)
    a = 0;
end
if nargin < 4 || isempty(b)
    b = 0;
end

% check inputs
if any( a(:) < 0 ) || ~isscalar(a) && any(size(a) ~= size(x))
    error('''a'' must be non-negative and either scalar or the same size as ''x''');
end
if ~isscalar(b) && any(size(b) ~= size(x))
    error('''b'' must be scalar or the same size as ''x''');
end
%-----%


% compute the prox
y = prox_phi( (x - gamma .* b) ./ (2*a*gamma + 1), gamma ./ (2*a*gamma + 1), varargin{:} );