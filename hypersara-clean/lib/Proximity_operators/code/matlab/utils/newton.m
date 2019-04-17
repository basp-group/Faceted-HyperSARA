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


function p = newton(fun, der, p, low, high)

% default parameters
TOL = 1e-7;
MAX = 100;
if nargin < 4, low  = eps; end
if nargin < 5, high = Inf; end
%-----%

for it = 1:MAX
    
    p_old = p;
    
    p = p - fun(p) ./ der(p);  % newton step
    p = min(max(low,p),high);  % range constraint
    
    err = abs(p-p_old) ./ abs(p_old);
    
    if all(err(:) <= TOL)
        break;
    end
    
end