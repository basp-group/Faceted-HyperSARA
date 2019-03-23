 function p = project_L12(x, eta, dir)
%function p = project_L12(x, eta, dir)
%
% This procedure computes the projection onto the constraint set:
%
%              ||x^(1)||_2 + ... + ||x^(B)||_2 <= eta
%
%  INPUTS
% ========
%  x   - ND array
%  eta - positive, scalar
%  dir - integer, direction along which the blocks are stored
% 
%  DEPENDENCIES
% ==============
%  project_L1.m - located in the folder 'indicator'
%  prox_Linf.m  - located in the folder 'multi'
%  prox_max.m   - located in the folder 'multi'

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

% compute the projection
if sum( xa(:) ) <= eta
    p = x;
else
           
    % project the block norms
    ya = project_L1(xa, eta);
    
    % compute the scale factors
    ya = ya ./ xa;
    ya(xa < eps) = 0;
    
    % rescace the blocks
    p = bsxfun(@times, x, ya);
    
end