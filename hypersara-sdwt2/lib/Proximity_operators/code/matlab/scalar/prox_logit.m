 function p = prox_logit(x, gamma)
%function p = prox_logit(x, gamma)
%
% This procedure computes the proximity operator of the function:
%
%            f(x) = gamma * log( 1 + exp(x) )
%
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'
% 
%  DEPENDENCIES
% ==============
%  prox_entropy_symm.m - located in the folder 'scalar'

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
%-----%


%Use Moreau's decomposition formula
p = x - gamma.*prox_entropy_symm(x./gamma,1./gamma);

% 
% limit = 5e2;
% 
% %RUN HALEY'S METHOD FOR SOLVING w = W_{exp(x)}(gamma*exp(x))
% %W_r, r real, is the generalized Lambert function defined in 
% %Mezo et al. "On the generalization of the Lambert W function", 
% %Tech. Rep., 2015, http://arxiv.org/abs/1408.3999
% 
% w = zeros(size(x));
% ex = exp(x);
% z = gamma.*ex;
%  
% %INITIALIZATION
% approx = gamma.*(1 - exp(gamma-x));
% w(z>1) = approx(z>1);
% 
% %RUN
% maxiter = 20;
% testend = zeros(size(x));
% prec = 1e-8;
% 
% for it = 1:maxiter
%  
%     e = exp(w);
%     y = w.*e + ex.*w - z; 
%     v = e.*(1 + w) + ex; 
%     u = e.*(2 + w);
%     wnew = w - y./(v - y.*u./(2*v));
%     testend(abs(wnew-w)./abs(w)  < prec & testend == 0) = 1;
%     idx_update = find( abs(wnew-w)./abs(w)  >= prec & testend == 0);
%     w(idx_update) = wnew(idx_update); %the rest stays constant !
%     if(sum(testend)==length(w))%stop !
%         break;
%     end
% 
% end
%  
% p = x - w;
% 
% %ASYMPTOTIC DVP
% p(x>limit) = x(x>limit) - approx(x>limit);
 
