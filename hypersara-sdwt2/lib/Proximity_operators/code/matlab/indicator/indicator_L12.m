 function p = indicator_L12(x, eta, dir)
%function p = indicator_L12(x, eta, dir)
%
% This procedure evaluates the indicator function of the constraint set:
%
%              ||x^(1)||_2 + ... + ||x^(B)||_2 <= eta
%
%  INPUTS
% ========
%  x   - ND array
%  eta - positive, scalar
%  dir - integer, direction along which the blocks are stored

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


% check input
if nargin < 2 || ~isscalar(eta) || eta <= 0
    error('''eta'' must be positive and scalar');
end
if nargin < 3 || isempty(dir) || ~isscalar(dir) || dir <= 0
    error('''dir'' must be a positive integer');
end
%-----%


% compute the block norms
xa  = sqrt( sum(x.^2,dir) );

% evaluate the indicator function
if sum( xa(:) ) <= eta
    p = 0;
else  
    p = Inf;
end