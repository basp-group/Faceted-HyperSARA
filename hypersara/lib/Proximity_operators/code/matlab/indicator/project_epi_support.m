 function [p,t] = project_epi_support(y, xi, a, b)
%function [p,t] = project_epi_support(y, xi, a, b)
%
% This procedure computes the projection onto the epigraph of 
%
%                        phi(y) = sigma_[a,b](y)
%
% When the inputs are arrays, the outputs are computed element-wise.
%
%  INPUTS
% ========
%  y  - ND array
%  xi - ND array with the same size as 'y'
%  a  - negative, scalar or ND array
%  b  - positive, scalar or ND array

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
if any(a(:) >= 0) || any(b(:) <= 0)
    error('''a'' must be negative and ''b'' must be positive');
end
if ~isscalar(a) && any(size(a) ~= size(y))
    error('''a'' must be either scalar or the same size as ''y''')
end
if ~isscalar(b) && any(size(b) ~= size(y))
    error('''b'' must be either scalar or the same size as ''y''')
end
%-----%


% 4th branch
p = zeros(size(y));
t = p;

% 3rd branch
mask = a.*y > xi & -y./a <= xi;
pp = (y+a.*xi) ./ (1+a.^2);
tt = a .* pp;
p(mask) = pp(mask);
t(mask) = tt(mask);

% 2nd branch
mask = b.*y > xi & -y./b <= xi;
pp = (y+b.*xi) ./ (1+b.^2);
tt = b .* pp;
p(mask) = pp(mask);
t(mask) = tt(mask);

% 1st branch
mask = a.*y <= xi & b.*y <= xi;
p(mask) =  y(mask);
t(mask) = xi(mask);