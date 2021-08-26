clc; clear all; close all;
format compact;

addpath ../../../CubeHelix
addpath ../../lib/print
addpath ../../results/mnras_faceted_corrected/final_real_data

%% Parameters
load_images = 1;
load_residuals = 1;
fig_size = [1250, 1250]; % only change wrt Abdullah's code
pixel_shift_colorbar = 30;
extension = '.pdf';
map_img = cubehelix(2048);
fontsize=60;

%% Load images
if load_images
    xhs = fitsread('x_fhs_reduced.fits');
    xl1 = fitsread('xl1_reduced.fits');
    xclean = fitsread('xclean_reduced.fits');
end
[N1, N2, c] = size(xhs);

xclean(:,:,1) = xclean(:,:,1) * 42.1827;
xclean(:,:,2) = xclean(:,:,2) * 8.3285;
xclean(xclean < 0) = 0;

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
psf_flux = [40.7 7.97 1];


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
%     % revise scale zoom (ticks are weird)
%     [f, h] = display_zoom(xclean_z1(:,:,band), fig_size, pixel_shift_colorbar, psf_flux(band).*clim_log, map_img, fontsize, round(psf_flux(band).*tickLabel,4)); % round(psf_flux(band).*[0.01,0.02,0.03,0.04],2)
%     export_fig(['figs/xclean_z11_ch', num2str(band),extension], '-transparent','-q101')
%     close
% end

%% Plot the west hot spot
clim_log = [0.0005,0.2]; % HS and SARA
tickLabel = [0.001,0.01,0.1];
%clim_log_cl = psf_flux * [0.008 0.05]; % CLEAN
%tickLabel_cl = psf_flux(band) * [0.01,0.04];
cen = [1110 422];
len = 120;

xhs_z1 = xhs(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
xl1_z1 = xl1(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
xclean_z1 = xclean(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);

for band = 1:size(xhs,3)
    [f, h] = display_zoom(xhs_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
    exportgraphics(f,['figs/xhs_z31_ch', num2str(band),extension], ...
    'ContentType', 'vector', 'BackgroundColor','none')
    close
    
    %
    [f, h] = display_zoom(xl1_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
    exportgraphics(f,['figs/xl1_z31_ch', num2str(band),extension], ...
    'ContentType', 'vector', 'BackgroundColor','none')
    close
     
    %
    [f, h] = display_zoom(xclean_z1(:,:,band), fig_size, pixel_shift_colorbar, psf_flux(band).*clim_log, map_img, fontsize, round(psf_flux(band).*tickLabel,4)); % round(psf_flux(band).*[0.01,0.02,0.03,0.04],2)
    exportgraphics(f,['figs/xclean_z31_ch', num2str(band),extension], ...
    'ContentType', 'vector', 'BackgroundColor','none')
    close
end

%% Plot the inner core
clim_log = [1e-5  0.02];
tickLabel = [1e-4,1e-3,1e-2]; % [1e-4,1e-3]; %[1e-4,1e-3,1e-2];
%clim_log_cl = psf_flux * [0.00005  0.0015]; % clean
%tickLabel_cl = round(psf_flux(band) .*[1e-4,1e-3]);
cen = [764 1298];
len = 57;

xhs_z1 = xhs(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
xl1_z1 = xl1(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);
xclean_z1 = xclean(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2,:);

for band = 1 : size(xhs,3)
    [f, h] = display_zoom(xhs_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
    exportgraphics(f,['figs/xhs_z21_ch', num2str(band),extension], ...
    'ContentType', 'vector', 'BackgroundColor','none')
    close
    
    %
    [f, h] = display_zoom(xl1_z1(:,:,band), fig_size, pixel_shift_colorbar, clim_log, map_img, fontsize, tickLabel);
    exportgraphics(f,['figs/xl1_z21_ch', num2str(band),extension], ...
    'ContentType', 'vector', 'BackgroundColor','none')
    close
    
    %
    [f, h] = display_zoom(xclean_z1(:,:,band), fig_size, pixel_shift_colorbar, psf_flux(band).*clim_log, map_img, fontsize, round(psf_flux(band).*tickLabel,4));
    % if band == 1
    %     set(h,'XTick',[3e-3,5e-2],'XTickLabel',{'3e-3','5e-2'});
    % elseif band == 2
    %     set(h,'XTick',[1e-3,1e-2],'XTickLabel',{'1e-3','1e-2'});
    % else
    %     set(h,'XTick',[1e-4,1e-3],'XTickLabel',{'1e-4','1e-3'});
    % end
    exportgraphics(f,['figs/xclean_z21_ch', num2str(band),extension], ...
    'ContentType', 'vector', 'BackgroundColor','none')
    close
end