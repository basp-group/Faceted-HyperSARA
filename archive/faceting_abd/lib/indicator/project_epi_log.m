 function [p,t] = project_epi_log(y, xi, w)
%function [p,t] = project_epi_log(y, xi, w)
%
% This procedure computes the projection onto the epigraph of 
%
%                      phi(y) = -w * log(y)
%
% When the inputs are arrays, the outputs are computed element-wise.
%
%  INPUTS
% ========
%  y  - ND array
%  xi - ND array with the same size as 'y'
%  w  - positive, scalar or ND array

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


% check inputs
if nargin < 3 || isempty(w)
    w = 1;
end
if any(w(:) <= 0) || ~isscalar(w) && any(size(w) ~= size(y))
    error('''w'' must be positive and either scalar or the same size as ''y''')
end
%-----%

TOL = 1e-7;
MAX = 100;

lambda = eps;
for it = 1:MAX
	
	p = 0.5 * ( y + sqrt(y.^2 + 4*w.*lambda) );
	t = xi + lambda;
	
	d = -w.*log(p+eps) - t;
	h = -w.^2 ./ (p.^2 + w.*lambda) - 1;
    
    la_old = lambda;
	lambda = max(0, lambda - d ./ h);
	
    err = abs(lambda - la_old) ./ abs(la_old);
    
    if all(err(:) <= TOL)
        break;
    end
end