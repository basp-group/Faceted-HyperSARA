 function y = prox_svd(x, gamma, prox_phi, varargin)
%function y = prox_svd(x, gamma, prox_phi, varargin)
%
% This procedure computes the proximity operator of the matrix function:
%
%                    f(X) = gamma * phi(s)     with   X = U' * diag(s) * V
% 
% It is assumed that the matrices are stored along the dimensions 1 and 2.
%
%  INPUTS
% ========
%  x       - ND array
%  gamma   - positive, scalar or ND array compatible with the size of 'x'
% prox_phi - function handle with two arguments at least
% varargin - additional parameters for the function 'prox_phi' [OPTIONAL]
% 
%  DEPENDENCIES
% ==============
%  sv_dec.cpp     - [MEX-FILE] located in the folder 'utils'
%  sv_dec_fat.cpp - [MEX-FILE] located in the folder 'utils'
%  sv_rec.cpp     - [MEX-FILE] located in the folder 'utils'

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

% prox computation
g = prox_phi(s, gamma, varargin{:});

% spectral reconstruction
y = sv_rec(u, g, v);