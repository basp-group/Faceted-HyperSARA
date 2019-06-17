 function y = fun_svd(x, gamma, fun_phi, varargin)
%function y = fun_svd(x, gamma, fun_phi, varargin)
%
% This procedure evaluates the function:
%
%                f(X) = gamma * phi(s)     with   X = U' * diag(s) * V
% 
% It is assumed that the matrices are stored along the dimensions 1 and 2.
%
%  INPUTS
% ========
%  x       - ND array
%  gamma   - positive, scalar or ND array compatible with the size of 'x'
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


% spectral decomposition
if size(x,2) > size(x,1)
    [u,s,v] = sv_dec_fat(x);
else
    [u,s,v] = sv_dec(x);
end

% function computation
y = fun_phi(s, gamma, varargin{:});