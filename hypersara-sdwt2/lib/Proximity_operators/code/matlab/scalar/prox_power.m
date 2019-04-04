 function p = prox_power(x, gamma, q)
%function p = prox_power(x, gamma, q)
%
% This procedure computes the proximity operator of the function:
%
%                   f(x) = gamma * |x|^q     with q >= 1
%
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
%  q     - not less than 1, scalar
% 
%  DEPENDENCIES
% ==============
%  newton.m - located in the folder 'utils'
%
%  REMARKS
% =========
% For many values of q, the computation is based on Newton's method.
% Depending on the inputs, you may experience slowness or numerical issues.

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
if any( q(:) < 1 ) || ~isscalar(q)
    error('''q'' must be a real number not less than 1')
end
%-----%


% compute the prox
if q == 1
    
    p = sign(x) .* max(0, abs(x) - gamma);
    
elseif q == 4/3
    
    xi = sqrt(x.^2 + 256 * gamma.^3 / 729);
    p = x + 4*gamma .* (nthroot(xi-x,3) - nthroot(xi+x,3)) / (3*nthroot(2,3));
	
elseif q == 3/2

    p = x + 9/8 * gamma.^2 .* sign(x) .* ( 1 - sqrt(1 + 16/9 * abs(x) ./ gamma.^2) );
					
elseif q == 2
    
    p = x ./ (2*gamma+1);
    
elseif q == 3
    
	p = sign(x) .* ( sqrt(1 + 12*gamma .* abs(x)) - 1 ) ./ (6*gamma);
                    
elseif q == 4

    xi = sqrt( x.^2 + 1./(27*gamma) );
    p = nthroot((xi+x)./(8*gamma), 3) - nthroot((xi-x)./(8*gamma), 3);
                 
else
    
    % prepare the Newton's method
    fun = @(t) q * gamma .* t.^(q-1) + t - abs(x);
    der = @(t) q*(q-1) * gamma .* t.^(q-2) + 1;
    
    % initialize the solution
    p_init = abs(x);
    
    % use the Newton's method
    p = newton(fun, der, p_init);
    
    % compute the prox   
    p = sign(x) .* p;
                    
end