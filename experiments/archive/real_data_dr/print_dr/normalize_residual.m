clc; clear all; close all;

addpath ../../mnras_faceted_corrected/final_dr/gamma5e-6;
addpath ../../mnras_faceted_corrected/final_dr/;
load('psf_normalization.mat');
% psf_flux = [flux(1:2), 1];

%%
% % renormalize residual images (DR)
% res = fitsread('facethyperSARA_res_it7000_gamma5e-06_gamma0_0.0001_2b_fouRed2_perc15_adpteps0.fits');
% res = res./reshape(scale_factors, [1, 1, 30]);
% fitswrite(res, 'res_hs_1e-4.fits');
%
% %%
% res = fitsread('facethyperSARA_res_it7000_gamma5e-06_gamma0_0.0005_2b_fouRed2_perc15_adpteps0.fits');
% res = res./reshape(scale_factors, [1, 1, 30]);
% fitswrite(res, 'res_hs_5e-4.fits');

%%
% res = fitsread('facethyperSARA_res_it10752_reweight20_gamma5e-06_gamma0_0.001_2b_fouRed2_perc15_adpteps0.fits');
res = fitsread('facethyperSARA_res_it10500_reweight20_gamma5e-06_gamma0_0.0005_2b_fouRed2_perc15_adpteps0.fits');
res = res ./ reshape(scale_factors, [1, 1, 30]);
r = zeros(size(res, 1), size(res, 2), 3);
r(:, :, [1, 2]) = res(:, :, [1, end]);
r(:, :, 3) = mean(res, 3);
fitswrite(res, 'r_fhs_dr_5e-4.fits');

%% Reading only a portion of a fits file

% extract first and last channels from the fits file
% info        = fitsinfo('facethyperSARA_res_it7000_gamma5e-06_gamma0_0.001_2b_fouRed2_perc15_adpteps0.fits');
% rowend      = info.PrimaryData.Size(1);
% colend      = info.PrimaryData.Size(2);
% sliceend    = info.PrimaryData.Size(3);
% res = fitsread('facethyperSARA_res_it7000_gamma5e-06_gamma0_0.001_2b_fouRed2_perc15_adpteps0.fits','primary',...
%           'Info', info,...
%           'PixelRegion',{[1 1 rowend], [1 1 colend], [1 sliceend-1 sliceend]});
