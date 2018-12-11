 function p = prox_inner(x, gamma, w, z, prox_phi, varargin)
%function p = prox_inner(x, gamma, w, z, prox_phi, varargin)
%
% This procedure computes the proximity operator of a function:
% 
%                  f(x) = gamma * phi(x/w - z)
%
% When the input 'x' is an array, the expression 'x/w - z' is meant element-wise.
%
%  INPUTS
% ========
% x        - ND array
% gamma    - positive, scalar or ND array with the same size as 'x'
% w        - non-zero, scalar or ND array with the same size as 'x' [OPTIONAL]
% z        - scalar or ND array with the same size as 'x' [OPTIONAL]
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
if nargin < 3 || isempty(w)
    w = 1;
end
if nargin < 4 || isempty(z)
    z = 0;
end

% check inputs
if any( w(:) == 0 ) || ~isscalar(w) && any(size(w) ~= size(x))
    error('''w'' must be non-zero and either scalar or the same size as ''x''');
end
if ~isscalar(z) && any(size(z) ~= size(x))
    error('''z'' must be scalar or the same size as ''x''');
end
%-----%


% compute the prox
p = prox_phi(x./w - z, gamma ./ w.^2, varargin{:});
p = w .* (z + p);