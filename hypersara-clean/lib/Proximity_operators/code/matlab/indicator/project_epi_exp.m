 function [p,t] = project_epi_exp(y, xi, w)
%function [p,t] = project_epi_exp(y, xi, w)
%
% This procedure computes the projection onto the epigraph of 
%
%                      phi(y) = w * exp(y)
%
% When the inputs are arrays, the outputs are computed element-wise.
%
%  INPUTS
% ========
%  y  - ND array
%  xi - ND array with the same size as 'y'
%  w  - positive, scalar or ND array
% 
%  DEPENDENCIES
% ==============
%  newton.m - located in the folder 'utils'

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
    error('''gamma'' must be positive and either scalar or the same size as ''y''')
end
%-----%


% 1st branch
p = y;
t = max(xi, w .* exp(p));

% branch selection
mask = (w .* exp(y) > xi) & (y > log(eps)+1);
yy = y(mask);
xx = xi(mask);
if isscalar(w)
    ww = w;
else
    ww = w(mask);
end
    
% newton's method
fun = @(t)   t.^2 - xx .* t + log(t) - yy;
der = @(t) 2*t    - xx      + 1./t;

% root finding
rho = max(xx,eps);
rho = newton(fun, der, rho);

% 2nd branch
p(mask) = log(rho./ww);
t(mask) = rho;