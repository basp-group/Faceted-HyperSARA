function [p,q] = prox_conicl0(x,d,gamma,Delta)
%function [p,q] = prox_conicl0(x,d,gamma,Delta)
%
% This procedure computes the proximity operator of the function:
%
%                f(x,d) = gamma  \ell_0(x) + \iota_{S}(x,d) 
%
% with S = {x \in C^N, d \in C^N, s.t. (\forall n \in {1,...,N}) 
%           \exist delta_n \in [-Delta_n,Delta_n] s.t. d_n = \delta_n x_n 
%           }
% where \Delta \in [0,+\infty)^N
%
% When the inputs '(x,d)' are arrays, the outputs '(p,q)' are computed element-wise.
%
%  INPUTS
% ========
%  x     - ND array complex valued
%  d     - ND array complex valued with the same size as 'x'
%  gamma - positive, scalar 
%  Delta - positive or null real valued ND array with the same size as 'x'

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
if ((~isscalar(gamma))||(gamma<=0))
    error('''gamma'' must be positive, scalar or the same size as ''x''')
end
if any( Delta(:) < 0 ) || any(size(Delta) ~= size(x))
    error('''Delta'' must be positive and with the same size as ''x''')
end
if any(size(d) ~= size(x))
    error('''d'' must be positive and with the same size as ''x''')
end
 
p = zeros(size(x));
q = zeros(size(x));

z = real(x.*conj(d));

eta = sqrt(((abs(d).^2 - abs(x).^2)).^2 + 4*z.^2);

deltahat = Delta.*ones(size(x));
deltahat(z ~=0) = min((eta + abs(d).^2 - abs(x).^2)./(2*abs(z)),Delta).*sign(z);
deltahat(z == 0 & abs(x) >= abs(d)) = 0;
 
seuilmod = (abs(deltahat.*x - d).^2)./(1 + deltahat.^2) + 2*gamma;
indseuil = find(abs(x).^2 + abs(d).^2 >= seuilmod);
 
upmod = (x + deltahat.*d)./(1 + deltahat.^2);

p(indseuil) = upmod(indseuil);
q(indseuil) = upmod(indseuil).*deltahat(indseuil);