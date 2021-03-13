clc; clear all; close all;
format compact;

addpath ../../../CubeHelix
addpath ../../../export_fig_master
addpath ../../lib/print
addpath ../../mnras_faceted_corrected/final_dr/

%% Parameters
load_images = 1;
load_residuals = 1;
fig_size = [1250, 1250];
shift_colorbar = [0,0,0,0]; % + [left, bottom, width, height] to place it where you want
pixel_shift_colorbar = 30;
extension = '.pdf';
map_img = cubehelix(2048);
fontsize=60;
load('flux_psf.mat')
psf_flux = [flux(1:2), 1];

%% Load images
if load_images
    xhs = fitsread('x_fhs_reduced_dr_5e-4.fits');
    xhs_avg = fitsread('x_fhs_avg_reduced.fits');
    xl1 = fitsread('xl1_reduced_dr.fits');
    xclean = fitsread('xclean_reduced_dr.fits');
end
[N1, N2, c] = size(xhs);
xhs = flipud(xhs);
xhs_avg = flipud(xhs_avg);

% mutliply by l1 norm of clean-beam, contained in flux
xclean(xclean < 0) = 0;
xclean = xclean.*reshape(psf_flux, [1, 1, 3]);

%% Take only the effective window of the image
% a1 = 200; %up
% a2 = 250; %down
% b1 = 220; %left
% b2 = 100;
% 
% xhs1 = xhs(a1:end-a2,b1:end-b2,:);
% xl11 = xl1(a1:end-a2,b1:end-b2,:);
% xclean1 = xclean(a1:end-a2,b1:end-b2,:);
% 
% rhs1 = rhs(a1:end-a2,b1:end-b2,:);
% rl11 = rl1(a1:end-a2,b1:end-b2,:);
% rclean1 = rclean(a1:end-a2,b1:end-b2,:);

%% Diplay reasults (images, just first and last band)
% ground truth, faceted hyperSARA, hyperSARA, SARA

%=========================================================================%
% Plot parameters
%=========================================================================%
mkdir figs

%% Plot the east hot spot
% clim_log = [0.0005,0.14]; % HS and SARA
% tickLabel = [0.001,0.01,0.1];
% %clim_log_cl = psf_flux * [0.008 0.05]; % CLEAN
% %tickLabel_cl = psf_flux(band) * [0.01,0.04];
% cen = [410 2282];
% len = 130;
% 
% xhs_z1 = xhs(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
% xl1_z1 = xl1(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
% xclean_z1 = xclean(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
% 
% for band = 1:size(xhs,3)
% 
%     [f, h] = display_zoom(xhs_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
%     export_fig(['figs/xhs_z11_ch', num2str(band),extension],'-transparent','-q101')
%     close
%     
%     %
%     [f, h] = display_zoom(xl1_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
%     export_fig(['figs/xl1_z11_ch', num2str(band),extension], '-transparent','-q101')
%     close
% 
%     %
%     % ...
% end

%% Plot the west hot spot
clim_log = [0.0005,0.2]; % HS and SARA
tickLabel = [0.001,0.01,0.1];
%clim_log_cl = psf_flux * [0.008 0.05]; % CLEAN
%tickLabel_cl = psf_flux(band) * [0.01,0.04];
cen = [1110 422];
len = 120;

xhs_z1 = xhs(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
xhs_avg_z1 = xhs_avg(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
xl1_z1 = xl1(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
xclean_z1 = xclean(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);

for band = 1:size(xhs,3)
    [f, h] = display_zoom(xhs_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
    pause(0.5)
    export_fig(['figs/xhs_z31_ch', num2str(band),extension],'-transparent','-q101')
    
    [f, h] = display_zoom(xhs_avg_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
    pause(0.5)
    export_fig(['figs/xhs_avg_z31_ch', num2str(band),extension],'-transparent','-q101')
    
    [f, h] = display_zoom(xl1_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
    pause(0.5)
    export_fig(['figs/xl1_z31_ch', num2str(band),extension], '-transparent','-q101')

    [f, h] = display_zoom(xclean_z1(:,:,band), fig_size, pixel_shift_colorbar, psf_flux(band).*clim_log, map_img, fontsize, round(psf_flux(band) .* tickLabel,2));
    pause(0.5)
    export_fig(['figs/xclean_z31_ch', num2str(band),extension], '-transparent','-q101')
end
close all

%% Plot the inner core
clim_log = [1e-5  0.02];
tickLabel = [1e-4,1e-3,1e-2];
%clim_log_cl = psf_flux * [0.00005  0.0015]; % clean
%tickLabel_cl = round(psf_flux(band) .*[1e-4,1e-3]);
cen = [764 1298];
len = 57;

xhs_z1 = xhs(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
xhs_avg_z1 = xhs_avg(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
xl1_z1 = xl1(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
xclean_z1 = xclean(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);

for band = 1 : size(xhs,3)
    [f, h] = display_zoom(xhs_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
    pause(0.5)
    export_fig(['figs/xhs_z21_ch', num2str(band),extension], '-transparent','-q101')
    
    [f, h] = display_zoom(xhs_avg_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
    pause(0.5)
    export_fig(['figs/xhs_avg_z21_ch', num2str(band),extension], '-transparent','-q101')
    
    [f, h] = display_zoom(xl1_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
    pause(0.5)
    export_fig(['figs/xl1_z21_ch', num2str(band),extension], '-transparent','-q101')

    [f, h] = display_zoom(xclean_z1(:,:,band), fig_size, pixel_shift_colorbar, psf_flux(band) .* clim_log, map_img, fontsize, round(psf_flux(band) .* tickLabel,4));
    pause(0.5)
    export_fig(['figs/xclean_z21_ch', num2str(band),extension], '-transparent','-q101')
end
close all
