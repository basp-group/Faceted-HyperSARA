function [xyz] = baseline2xyz(d, az, elev, lat)
%%code to simulate continuous uv coverage with tracking from
% http://www.astro.umd.edu/~cychen/MATLAB/ASTR410/uvAndBeams.html#5

% Function computes transformation from baseline to x, y, z coordinates
% Inputs:
%   d: baseline length
%   az: baseline azimuth angle (E from N) [radians]
%   elev: baseline elevation angle [radians]
%   lat: observatory latitude [radians]
% Output:
%   vector containing xyz vector, length units same as baseline length
% Reference: Thompson, Moran and Swenson, eq. 4.4
%
% AH 2010.3.16

xyz = d * [cos(lat)*sin(elev) - sin(lat)*cos(elev)*cos(az);
           cos(elev)*sin(az);
           sin(lat)*sin(elev) + cos(lat)*cos(elev)*cos(az)];
       