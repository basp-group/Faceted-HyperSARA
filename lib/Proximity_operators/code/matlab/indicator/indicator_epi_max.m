 function p = indicator_epi_max(y, xi, w, dir)
%function p = indicator_epi_max(y, xi, w, dir)
%
% The procedure evaluates the indicator function of the epigraph of
%
%                     phi(y) = w * max_i y_i
%
% When the input 'y' is an array, the computation can vary as follows:
%  - dir = 0 --> 'y' is processed as a single vector [DEFAULT]
%                (in this case, 'xi' must be scalar)
%  - dir > 0 --> 'y' is processed block-wise along the specified direction
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


% linearize
if isempty(dir)
    y = y(:);
    dir = 1;
end

% check the constraint
mask = w .* max(y, [], dir) <= xi;

% evaluate the indicator function
if all(mask(:))
	p = 0;
else
	p = Inf;
end