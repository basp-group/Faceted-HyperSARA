 function p = project_monotone(x, dir)
%function p = project_monotone(x, dir)
%
% This procedure computes the projection onto the constraint set:
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
%
%  DEPENDENCIES
% ==============
%  pava.cpp - [OPTIONAL, MEX-FILE] located in the folder 'utils'
%
%  REMARKS
% =========
% The mex-file 'pava.cpp' provides a C++ implementation of the *pool adjacency
% violator algorithm* (PAVA). When the latter is not present in MATLAB's path,
% this procedure executes a (SLOWER) matlab implementation.
%
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


% linearize
sz = size(x);
if isempty(dir)
    x = x(:);
end

% select the implementation
if exist('pava', 'file')
    x = reshape(x, [sz(1) prod(sz(2:end))]);
    p = pava(x, ones(size(x)));
    p = reshape(p, sz);
else
    p = pava_matlab(x);
end





function p = pava_matlab(x)

S = cumsum_matrix(x);
S = cummin(S,2,'reverse');
p = max(S,[],1);
p = reshape(p, size(x));

function S = cumsum_matrix(x)
%
% This function computes the sum:
%
%    S(j,k,:) = sum( x(j:k,:) ) / (k-j+1)
%
% for every 'j' and 'k' in {1,...,N}, where N = size(x,1).


sz = size(x);
S = -Inf([sz(1) sz]);

s = cumsum([zeros([1 sz(2:end)]); x]);
s = reshape(s, [1 sz(1)+1 sz(2:end)]);

for i=1:sz(1)
    S(i,i:end,:) = bsxfun( @minus, s(:,i+1:end,:), s(:,i,:) );
    S(i,i:end,:) = bsxfun( @rdivide, S(i,i:end,:), 1:sz(1)-i+1 );
end
