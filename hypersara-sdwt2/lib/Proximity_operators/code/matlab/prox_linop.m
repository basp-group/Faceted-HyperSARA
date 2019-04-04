 function p = prox_linop(x, gamma, A, b, prox_phi, varargin)
%function p = prox_linop(x, gamma, A, b, prox_phi, varargin)
%
% This procedure computes the proximity operator of a function:
%
%                   f(x) = gamma * phi(A x + b)
% 
% where A is a semi-orthogonal linear operator (i.e., AA' = nu I).
%
%  INPUTS
% ========
% x        - ND array
% gamma    - positive, scalar
% A        - struct with the fields 'dir_op', 'adj_op', and 'nu'.
% b        - scalar or ND array with the same size as 'x' [OPTIONAL]
% prox_phi - function handle with two arguments at least
% varargin - additional parameters for the function 'prox_phi' [OPTIONAL]

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
if nargin < 4 || isempty(b)
    b = 0;
end

% check inputs
if ~isscalar(gamma)
    error('''gamma'' must be scalar');
end
if ~isstruct(A) || ~isfield(A,'dir_op') || ~isfield(A,'adj_op') || ~isfield(A,'nu')
    error('''A'' must be a struct with the fields ''dir_op'', ''adj_op'', and ''nu''');
end
if ~isscalar(b) % && any( size(b) ~= size(A.dir_op(x)) )
    error('''b'' must be scalar or the same size as ''Ax''');
end
%-----%


% compute the prox
p = prox_phi(A.dir_op(x) + b, A.nu * gamma, varargin{:});
p = x + 1/A.nu * A.adj_op(p - A.dir_op(x) - b);