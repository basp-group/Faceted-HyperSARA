 function p = project_log(x, eta, z, alpha)
%function p = project_log(x, eta, z, alpha)
%
% This procedure computes the projection onto the constraint set:
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

% target function
KLD = alpha * x - z + z .* log(z ./ (alpha.*x));

% projection    
if sum(KLD(:)) > eta
    
    % convert the upper bound
   lambda = get_KLD_lambda( x(:), z(:), alpha, eta );
   
   % compute the prox
   a = x - alpha*lambda;
   p = 0.5 * (a + sqrt(a.^2 + 4*lambda*z));
   
else
    p = x;
end






function lambda = get_KLD_lambda(x, z, alpha, eta)


% prox of KLD (and its derivative)
 t = @(b) 0.5 * ((x - alpha*b) + sqrt((x - alpha*b).^2 + 4*b*z));
dt = @(b) (z - t(b)) ./ sqrt((x - alpha*b).^2 + 4*b*z);

% discrep. function (and its derivative)
 c = sum( z - z .* log(z) );
 f = @(b) sum( t(b) - z .* log(t(b)) ) - c - eta;
df = @(b) sum( dt(b) .* (1 - z ./ t(b)) );

% find the zero of f
max_it = 1000;
tol = 10^-7;
lambda = 1;
for n = 1:max_it
    lambda_old = lambda;
    
    % Newton method
    lambda = lambda - f(lambda)/df(lambda);
    
    if abs(lambda-lambda_old) < tol * abs(lambda_old)
        break;
    end
end
if n == max_it
    warning('Reached the max number of iterations');
end