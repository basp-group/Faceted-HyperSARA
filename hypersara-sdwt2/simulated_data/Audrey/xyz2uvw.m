function [A] = xyz2uvw(h, delta)
%%code to simulate continuous uv coverage with tracking from
% http://www.astro.umd.edu/~cychen/MATLAB/ASTR410/uvAndBeams.html#5

% Function computes transformation from x, y, z to u, v, w coordinates
% Inputs:
%   h: sounce hour angle [radians]
%   delta: source declination [radians]
% Output:
%   matrix computing transform so uvw = A*xyz
% Reference: Thompson, Moran and Swenson, eq. 4.1
%
% AH 2010.3.16

A = [sin(h),             cos(h),             0;
    -sin(delta)*cos(h),  sin(delta)*sin(h),  cos(delta);
     cos(delta)*cos(h), -cos(delta)*sin(h),  sin(delta)];