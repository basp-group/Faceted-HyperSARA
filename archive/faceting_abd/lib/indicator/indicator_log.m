 function p = indicator_log(x, eta, z, alpha)
%function p = indicator_log(x, eta, z, alpha)
%
% This procedure evaluates the indicator function of the constraint set:
%
%    \sum_{n=1}^N  alpha * x_n - z_n + z_n * log(z_n/(alpha*x_n)) <= eta 
%
% Note that the input 'x' is processed as a single vector.
%
%  INPUTS
% ========
% x     - positive, ND array
% eta   - scalar
% z     - positive, ND array with the same size as 'x'
% alpha - positive, scalar [OPTIONAL]

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


% defaut inputs
if nargin < 4
    alpha = 1;
end

% check inputs
if any(x(:)) <= 0 || any(z(:)) <= 0 || any(alpha(:)) <= 0
    error('The inputs must be positive');
end
%-----%

% compute the KL divergence
KLD = alpha * x - z + z .* log(z ./ (alpha.*x));

% evaluate the indicator function    
if sum(KLD(:)) <= eta
    p = 0;
else
    p = Inf;
end