 function p = prox_root(x, gamma, q)
%function p = prox_root(x, gamma, q)
%
% This procedure computes the proximity operator of the function:
%
%                f(x) = gamma * |x|^q    with 0 < q < 1
%
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  q     - 0 <*< 1, scalar or ND array with the same size as 'x'
% 
%  DEPENDENCIES
% ==============
%  newton.m - located in the folder 'utils'
%
%  REMARKS
% =========
% Except for q=1/2, the computation is based on Newton's method. Hence,
% depending on the inputs, you may experience slowness or numerical issues.

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


% check input
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x''')
end
if any( q(:) <= 0 ) || any( q(:) >= 1 ) || ~isscalar(q) && any(size(q) ~= size(x))
    error('''q'' must belong to ]0,1[ and be either scalar or the same size as ''x''')
end
%-----%

% select the branches
alpha = gamma .* abs(x).^(q-2);
mask = alpha <= ( 2*(1-q)./(2-q) ).^(1-q) ./ (2-q);
alpha = alpha(mask);

% 1st branch
p = zeros(size(x));

% 2nd branch
if q == 0.5
   
    % explicit solution
    s = 2*sin( (acos(3*sqrt(3)*alpha/4) + pi/2) / 3 ) / sqrt(3);
    t = s.^2;
    
else
    
    % prepare the Newton's method
    fun = @(t) t - 1 + alpha .* q .* t.^(q-1);
    der = @(t) 1     + alpha .* q .* (q-1) .* t.^(q-2);
    
    % use the Newton's method
    t = newton(fun, der, 1);
    
end

% compute the prox
p(mask) = t .* x(mask);