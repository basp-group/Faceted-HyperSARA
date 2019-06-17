 function [p,t] = project_epi_distance(y, xi, w, dir, q, project_C, varargin)
%function [p t] = project_epi_distance(y, xi, w, dir, q, project_C, varargin)
%
% The procedure computes the projection onto the epigraph of
%
%                 phi(y) = w * ||y - P_C(y)||_2^q       with q=1 or q=2
%
% When the input 'y' is an array, the computation can vary as follows:
%
%  - dir = 0 --> 'y' is processed as a single vector [DEFAULT]
%                (in this case, 'xi' must be scalar)
%  - dir > 0 --> 'y' is processed block-wise along the specified direction
%                (in this case, 'xi' must be singleton along 'dir')
%  - dir < 0 --> 'y' is processed element-by-element.
%                (in this case, 'xi' must be the same size as 'y')
%
% The caller must ensure that 'project_C' uses the same type of processing.
%
%  INPUTS
% ========
%  y         - ND array
%  xi        - ND array compatible with the blocks of 'y'
%  w         - positive, scalar or ND array with the same size as 'xi'
%  dir       - integer, direction of block-wise processing [DEFAULT: 0]
%  q         - integer, scalar [DEFAULT: 2]
%  project_C - function handle with an argument at least [DEFAULT: @(y) 0]
%  varargin  - additional parameters for the function 'prox_phi' [OPTIONAL]
%
%  DEPENDENCIES
% ==============
%  solver3.m - located in the folder 'utils'
%
%  REMARKS
% =========
% The default values of 'q' and 'project_C' yield the projection onto the 
% epigraph of phi(y) = w * ||y||_2^2.

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
if nargin < 5 || isempty(q)
    q = 2;
end
if nargin < 6 || isempty(project_C)
    project_C = @(y) 0;
end

% check inputs
sz = size(y); if dir>0, sz(dir) = 1; end
if ~isempty(dir) && any( sz ~= size(xi) ) || isempty(dir) && numel(xi) ~= 1
    error('The input ''xi'' is not compatible with the size of ''y''');
end
if any( w(:) <= 0 ) || ~isscalar(w) && any(size(w) ~= sz)
    error('''w'' must be positive and either scalar or the same size as ''xi''');
end
%-----%


% linearize
sz = size(y);
if isempty(dir)
    y = y(:);
elseif dir<0
     y = reshape(y, [1 sz]);
    xi = reshape(xi, [1 sz]);
    dir = 1;
end

% preliminary
py = project_C(y);
dy = sqrt( sum((y-py).^2,dir) );

% scaling
if q == 1
    t = max(0, dy+w.*xi) ./ (1+w.^2);
elseif q == 2
    t = solver3(2*w.^2, 0, 1-2*w.*xi, -dy);
else
    error('Only ''q=1'' and ''q=2'' are supported');
end

% 3rd branch
p = py + bsxfun(@times, t./dy, y-py);
t = w .* sqrt( sum((p-project_C(p)).^2,dir) ).^q;

% 2nd branch
mask = dy == 0 & xi < 0;
p(mask) = y(mask);
t(mask) = 0;

% 1st branch
mask = w .* dy.^q <= xi;
p(mask) =  y(mask);
t(mask) = xi(mask);

% revert back
p = reshape(p, sz);
if numel(y) == numel(xi)
    t = reshape(t, sz);
end