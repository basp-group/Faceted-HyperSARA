clc; clear all; close all;

addpath data

%%
x0 = fitsread('W28_512.fits');

x = zeros(size(x0));
x(x0 > 10^(-3.9)) = x0(x0 > 10^(-3.9));
imagesc(log10(x)); colorbar

fitswrite(x, 'data/W28_512_m39.fits')