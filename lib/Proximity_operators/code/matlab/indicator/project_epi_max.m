 function [p,t] = project_epi_max(y, xi, w, dir)
%function [p,t] = project_epi_max(y, xi, w, dir)
%
% The procedure computes the projection onto the epigraph of
%
%                     phi(y) = w * max_i y_i
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
if any( w(:) <= 0 ) || ~isscalar(w) && any(size(w) ~= sz)
    error('''w'' must be positive and either scalar or the same size as ''xi''');
end
%-----%


% sort the inputs
v = bsxfun(@times, y, w);
v = sort(v, dir);

% compute the comulative sums
s = cat( dir, w.^2 .* xi, flip(v,dir) );
s = cumsum(s, dir);
s = flip(s, dir);

% normalize the cumulative sums
s = s ./ bsxfun( @minus, w.^2 + size(v,dir) + 1, cumsum(ones(size(s)),dir) );

% find the "right" value of t
infty  = inf(size(xi));
v_low  = cat(dir, -infty, v);
v_high = cat(dir, v, +infty);
mask = v_low < s & s <= v_high;

% extract the "right" value of t
mask = mask .* reshape( 1:numel(mask), size(mask) );
n = max(mask, [], dir);
t = reshape( s(n), size(xi) );

% compute p
p = bsxfun(@min, y, t./w);