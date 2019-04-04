 function p = indicator_monotone(x, dir)
%function p = indicator_monotone(x, eta, dir)
%
% This procedure evaluates the indicator function of the constraint set:
%
%                  x(1) <= x(2) <= ... <= x(N)
%
% When the input 'x' is an array, the computation can vary as follows:
%  - dir = 0 --> 'x' is processed as a single vector [DEFAULT]
%  - dir > 0 --> 'x' is processed block-wise along the FIRST direction
%
%  INPUTS
% ========
%  x   - ND array
%  dir - integer, direction of block-wise processing (FIRST ONLY)

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


% default inputs
if nargin < 3 || (~isempty(dir) && dir == 0)
    dir = [];
end
%-----%


% linearize
if isempty(dir)
    x = x(:);
end
dir = 1;

% check the constraint
mask = all( x(1:end-1,:) <= x(2:end,:), dir);

% evaluate the indicator function
if all(mask(:))
	p = 0;
else
	p = Inf;
end