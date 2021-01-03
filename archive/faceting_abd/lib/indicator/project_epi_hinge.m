 function [p,t] = project_epi_hinge(y, xi, w)
%function [p,t] = project_epi_hinge(y, xi, w)
%
% This procedure computes the projection onto the epigraph of 
%
%                    phi(y) = w * max(y,0)
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
    error('''gamma'' must be positive and either scalar or the same size as ''y''')
end
%-----%


% 3rd branch
c = max(0, y + w .* xi) ./ (1 + w.^2);
p = c;
t = c .* w;

% 1st and 2nd branches
mask = (w .* max(y,0) <= xi) | (y <= 0);
p(mask) =  y(mask);
t(mask) = max( xi(mask), 0 );