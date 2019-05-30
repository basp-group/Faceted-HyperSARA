clc; clear all; close all
format compact;

addpath data

x = fitsread('W28_1024.fits');
N = size(x);
facet_size = 256;

% generate hanning tapering function
d = 10; % number of facets modified on the borders
w_tmp = hann(2*d);
wx = [w_tmp(1:d); ones(N(2)-2*d, 1); w_tmp(d+1:end)];
wy = [w_tmp(1:d); ones(N(1)-2*d, 1); w_tmp(d+1:end)];
w = wy.*(wx.');

% apply tapering function on the borders
x_temp = w.*x;
figure; imagesc(log10(x_temp), [-5, 0]); colorbar; axis equal

% extend with extra zeros 
x = zeros(N+facet_size); % half facet size added to both dims for 50% overlap
x(facet_size/2+1:end-facet_size/2, facet_size/2+1:end-facet_size/2) = x_temp;

figure; imagesc(log10(x), [-5, 0]); colorbar; axis equal

% save result
% fitswrite('padded_W28_1024.fits', x)
