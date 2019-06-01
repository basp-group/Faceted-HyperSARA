clc; clear all; close all
format compact;

addpath data

x = fitsread('W28_1024.fits');
N_ref = size(x);
p = 0.5; % percentage of overlap (w.r.t the size of each facet)
Qy = 3; % number of facets
Qx = 3;
Q = Qx*Qy;
N = round(N_ref./(1 - 2*p./[Qy, Qx]));
padding = floor(p.*N./[Qy, Qx]);

% determine the size of the border:
% N_new - 2*p*N_new/Q = N, along each dimension

% generate hanning tapering function
d = 10; % number of facets modified on the borders
w_tmp = hann(2*d);
wx = [w_tmp(1:d); ones(N_ref(2)-2*d, 1); w_tmp(d+1:end)];
wy = [w_tmp(1:d); ones(N_ref(1)-2*d, 1); w_tmp(d+1:end)];
w = wy.*(wx.');

% apply tapering function on the borders
x_temp = w.*x;
figure; imagesc(log10(x_temp), [-5, 0]); colorbar; axis equal

% extend with extra zeros 
x = zeros(N); % half facet size added to both dims for 50% overlap
x(padding(1)+1:end-padding(1), padding(2)+1:end-padding(2)) = x_temp;

figure; imagesc(log10(x), [-5, 0]); colorbar; axis equal

% save result
fitswrite(x, 'data/padded_W28_1024.fits')
