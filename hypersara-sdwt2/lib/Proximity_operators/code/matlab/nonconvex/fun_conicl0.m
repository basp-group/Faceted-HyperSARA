 function p = fun_conicl0(x, d, Delta, gamma)
%function p = fun_conicl0(x, d, Delta, gamma)
%
% This procedure evaluates the function:
%
%                    f(x) = gamma * |x|^0 + \iota_{S}(x,d) 
%
% with S = {x \in C, d \in C, s.t.
%           \exist delta \in [-Delta,Delta] s.t. d = \delta x 
%           }
% where \Delta \in [0,+\infty)
%
% When the inputs '(x,d)' are arrays, the outputs '(p,q)' are the
% element-wise sum.
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

% check inputs
if any( gamma(:) <= 0 ) || ~isscalar(gamma) && any(size(gamma) ~= size(x))
    error('''gamma'' must be positive and either scalar or the same size as ''x''')
end
if any( Delta(:) < 0 ) || any(size(Delta) ~= size(x))
    error('''Delta'' must be positive and with the same size as ''x''')
end
if any(size(d) ~= size(x))
    error('''d'' must be positive and with the same size as ''x''')
end
%-----%

% evaluate the l0 function
p = sum( gamma(:) .* (x(:) ~= 0) );

% check the constraints

if(any(~isreal(x(:).*conj(d(:))))) %x and d does not have the same arg
    p = inf;
else
    %x and d have the same arg
    if(max(abs(x(d~=0))/abs(d(d~=0))) > Delta) %delta not in [-Delta;Delta]
        p = inf;
    end
end


