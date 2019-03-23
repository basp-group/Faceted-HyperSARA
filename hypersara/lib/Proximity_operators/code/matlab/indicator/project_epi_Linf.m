 function [p,t] = project_epi_Linf(y, xi, w, dir)
%function [p,t] = project_epi_Linf(y, xi, w, dir)
%
% The procedure computes the projection onto the epigraph of
%
%                     phi(y) = w * ||y||_inf
%
% When the input 'y' is an array, the computation can vary as follows:
%  - dir = 0 --> 'y' is processed as a single vector [DEFAULT]
%                (in this case, 'xi' must be scalar)
%  - dir > 0 --> 'y' is processed block-wise along the specified direction
%                (in this case, 'xi' must be singleton along 'dir')
%
%  INPUTS
% ========
%  y   - ND array
%  xi  - ND array compatible with the blocks of 'y'
%  w   - positive, scalar or ND array with the same size as 'xi'
%  dir - integer, direction of block-wise processing
% 
%  DEPENDENCIES
% ==============
%  project_epi_max.m - located in the folder 'indicator'

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
if nargin < 3 || isempty(w)
    w = 1;
end
if nargin < 4 || (~isempty(dir) && dir == 0)
    dir = [];
end
%-----%


% compute t
abs_y = abs(y);
[~,t] = project_epi_max(abs_y, xi, w, dir);
t = max(0,t);

% compute p
p = sign(y) .* bsxfun(@min, abs_y, t./w);