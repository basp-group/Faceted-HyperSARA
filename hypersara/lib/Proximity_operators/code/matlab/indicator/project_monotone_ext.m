 function p = project_monotone_ext(x, dir)
%function p = project_monotone_ext(x, dir)
%
% This procedure computes the projection onto the constraint set:
%
%               0 <= x(1) <= x(2) <= ... <= x(N)
%
% When the input 'x' is an array, the computation can vary as follows:
%  - dir = 0 --> 'x' is processed as a single vector [DEFAULT]
%  - dir > 0 --> 'x' is processed block-wise along the FIRST direction
%
%  INPUTS
% ========
%  x   - ND array
%  dir - integer, direction of block-wise processing (FIRST ONLY)
%
%  DEPENDENCIES
% ==============
%  project_monotone.m - located in the folder 'indicator'
%
%  REMARKS
% =========
% Unlike the other m-files, the block-wise processing can be only performed 
% along the 1st dimension of the input array.

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

p = project_monotone(x, dir);
p = max(p, 0);