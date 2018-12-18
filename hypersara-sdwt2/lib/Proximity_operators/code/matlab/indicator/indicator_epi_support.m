 function p = indicator_epi_support(y, xi, a, b)
%function p = indicator_epi_support(y, xi, a, b)
%
% This procedure evaluates the indicator functio of the epigraph of 
%
%                        phi(y) = sigma_[a,b](y)
%
% When the inputs are arrays, the output is the element-wise sum.
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


% check the constraint
mask = a.*y <= xi & b.*y <= xi;

% evaluate the indicator function
if all(mask(:))
	p = 0;
else
	p = Inf;
end