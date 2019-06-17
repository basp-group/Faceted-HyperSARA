 function p = prox_entropy_symm(x, gamma)
%function p = prox_entropy_symm(x, gamma)
%
% This procedure computes the proximity operator of the function:
%
%          f(x) = gamma * ( x * log(x) + (1-x) * log(1-x) )
%
% When the input 'x' is an array, the output 'p' is computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array
%  gamma - positive, scalar or ND array with the same size as 'x'

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

  
limit = 8;  
%RUN HALEY'S METHOD FOR SOLVING w = W_{exp(x/gamma)}(exp(x/gamma)/gamma)
%where W_r(x) is the generalized Lambert function

w = zeros(size(x));
igamma = 1./gamma;
c = x.*igamma - log(gamma);
z = exp(c);
r = exp(x.*igamma); 
 
% ASYMPTOTIC APPROX
approx = 1 - exp((1-x).*igamma);

%INITIALIZATION
% Case 1: gamma <= 1/30
w(z>1 & gamma <= 1/30) = c(z>1 & gamma <= 1/30) - log( c(z>1 & gamma <= 1/30) ) ;
% Case 2: gamma > 1/30
w(z>1 & gamma > 1/30) = igamma.*approx(z>1 & gamma > 1/30);
 
%RUN
maxiter = 20;
testend = zeros(size(x));
prec = 1e-8;

for it = 1:maxiter
 
    e = exp(w);
    y = w.*e + r.*w - z; 
    v = e.*(1 + w) + r; 
    u = e.*(2 + w);
    wnew = w - y./(v - y.*u./(2*v));
    testend(abs(wnew-w)./abs(w)  < prec & testend == 0) = 1;
    idx_update = find( abs(wnew-w)./abs(w)  >= prec & testend == 0);
    w(idx_update) = wnew(idx_update); %the rest stays constant !
    if(sum(testend)==length(w))%stop !
        break;
    end

end
 
p = gamma.*w;

%ASYMPTOTIC DVP
p(c>limit & gamma > 1) = approx(c>limit & gamma > 1);

%FINAL TRESHOLD TO AVOID NUMERICAL ISSUES FOR SMALL GAMMA
p = min(p,1);