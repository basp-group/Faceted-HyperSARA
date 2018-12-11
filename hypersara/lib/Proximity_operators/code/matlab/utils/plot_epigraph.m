 function plot_epigraph(range, name, varargin)
%function plot_epigraph(range, name, varargin)
%
% This function plots the projection onto the epigraph of a scalar function.
%
%  INPUTS
% =======
%    range - vector of points on which the epigraphical projection is computed
%     name - name of the epigraphical projection to be evaluated
% varargin - comma-separated list of additional parameters to the epigraphical projection

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

[y,xi] = meshgrid(range);

% evaluation
[p,t] = feval(['project_epi_' name], y, xi, varargin{:});

% visualization
plot(y, xi, '.b'); hold on;
plot(p, t, 'or');  hold off;