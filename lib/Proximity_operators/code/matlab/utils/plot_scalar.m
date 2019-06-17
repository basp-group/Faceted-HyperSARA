 function plot_scalar(range, name, varargin)
%function plot_scalar(range, name, varargin)
%
% This function plots the proximity operator of a scalar function.
%
%  INPUTS
% =======
%    range - vector of points on which the proximity operator is computed
%     name - name of the proximity operator to be evaluated
% varargin - comma-separated list of additional parameters to the proximity operator

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


% default input
if isempty(range)
    range = -1:0.1:1;
end

% prox. evaluation
prox = feval(['prox_' name], range, varargin{:});

% visualization
plot(range, prox, 'linewidth', 2, 'DisplayName', name);

% legend
legend('off');legend('show')
