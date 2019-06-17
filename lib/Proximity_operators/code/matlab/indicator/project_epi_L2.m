 function [p,t] = project_epi_L2(y, xi, w, dir)
%function [p t] = project_epi_L2(y, xi, w, dir)
%
% The procedure computes the projection onto the epigraph of
%
%                     phi(y) = w * ||y||_2
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

% check inputs
sz = size(y); sz(dir) = 1;
if ~isempty(dir) && any( sz ~= size(xi) ) || isempty(dir) && numel(xi) ~= 1
    error('The input ''xi'' is not compatible with the blocks of ''y''');
end
if any( w(:) < 0 ) || ~isscalar(w) && any(size(w) ~= sz)
    error('''w'' must be positive and either scalar or the same size as ''xi''');
end
%-----%


% linearize
sz = size(y);
if isempty(dir)
    y = y(:);
end

% compute the scaling
yy = sqrt( sum(y.^2, dir) );
a = max(0, 1 + w .* xi ./ yy) ./ (1+w.^2); % 3rd branch
a(w .* yy <= xi) = 1;                      % 2nd branch
a(yy==0) = 0;                              % 1st branch

% compute the projection
p = bsxfun(@times, a, y);
t = max(xi, w .* a .* yy);

% revert back
p = reshape(p, sz);
