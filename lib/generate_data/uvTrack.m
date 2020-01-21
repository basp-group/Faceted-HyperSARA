function [UVW] = uvTrack(h, d, az, el, lat, dec, nstep) 
%%code to simulate continuous uv coverage with tracking from
% http://www.astro.umd.edu/~cychen/MATLAB/ASTR410/uvAndBeams.html#5

%% Calculate a uv track for a single antenna pair
%
% Inputs
%   h;  hour angles [rad]
%   d;  baseline length
%   az; baseline azimuth angle [rad]
%   el; baseline elevation angle [rad]
%   lat; observatory latitude [rad]
%   dec; source declination [rad]
%   nstep; number of steps in calculation
% Output
%   matrix with u, v, w in columns. 
%
% AH 2010.3.16

%% Set up for calculation

% results matrices
UVW = zeros(nstep, 3); % u, v, w columns

%% Calculations
for i = 1:nstep
%     baseline2xyz(d, az, el, lat)
    UVW(i,:) = (xyz2uvw(h(i), dec) * baseline2xyz(d, az, el, lat))';
end

