clc ; clear all; close all;

addpath ../../mnras_faceted_corrected/final_dr/gamma5e-6
addpath ../../mnras_faceted_corrected/final_dr/
load('psf_normalization.mat')

%%

im = fitsread('facethyperSARA_xsol_it10500_reweight20_gamma5e-06_gamma0_0.0005_2b_fouRed2_perc15_adpteps0.fits');
m = mean(im, 3);
x = zeros(size(im, 1),size(im,2),3);
x(:,:,1:2) = im(:,:,[1, end]);
x(:,:,3) = m;
fitswrite(x, 'x_fhs_reduced_dr_5e-4.fits');

%%
res = fitsread('facethyperSARA_res_it10500_reweight20_gamma5e-06_gamma0_0.0005_2b_fouRed2_perc15_adpteps0.fits');
res = res./reshape(scale_factors, [1, 1, 30]);
r = zeros(size(res, 1), size(res, 2), 3);
r(:,:,[1,2]) = res(:,:,[1,end]);
r(:,:,3) = mean(res, 3);
fitswrite(res, 'r_fhs_reduced_dr_5e-4.fits');