 function p = project_hyperslab(x, low, high, w, dir)
%function p = project_hyperslab(x, low, high, w, dir)
%
% This procedure computes the projection onto the constraint set:
%
%                    low  <=  w'x  <=  high
% 
% When the input 'x' is an array, the computation can vary as follows:
%  - dir = 0 --> 'x' is processed as a single vector [DEFAULT]
%  - dir > 0 --> 'x' is processed block-wise along the specified direction
%
%  INPUTS
% ========
%  x    - ND array
%  low  - scalar or ND array compatible with the blocks of 'x' [OPTIONAL]
%  high - scalar or ND array compatible with the blocks of 'x' [OPTIONAL]
%  w    - ND array of the same size as 'x' [OPTIONAL]
%  dir  - integer, direction of block-wise processing [OPTIONAL]

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
if nargin < 2 || isempty(low)
    low = -inf;
end
if nargin < 3 || isempty(high)
    high = +inf;
end
if nargin < 4 || isempty(w)
    w = ones(size(x));
end
if nargin < 5 || (~isempty(dir) && dir == 0)
    dir = [];
end

% check inputs
sz = size(x); sz(dir) = 1;
if ~isscalar(low) && any(size(low) ~= sz)
    error('''low'' must be either scalar or compatible with the blocks of ''x''')
end
if ~isscalar(high) && any(size(high) ~= sz)
    error('''high'' must be either scalar or compatible with the blocks of ''x''')
end
if any( low(:) >= high(:) )
    error('''low'' must be lower than ''high''')
end
if any( size(w) ~= size(x) )
    error('''w'' must be the same size as ''x''')
end
%-----%


% linearize
sz = size(x);
if isempty(dir)
    x    = x(:);
    w    = w(:);
end

% prepare the offsets
xTw = sum(w .* x, dir);
wTw = sum(w .^ 2, dir);
mu1 = max(0, low  - xTw) ./ wTw;
mu2 = min(0, high - xTw) ./ wTw;

% compute the projection
p = x + bsxfun(@times, mu1 + mu2, w);

% revert back
p = reshape(p, sz);