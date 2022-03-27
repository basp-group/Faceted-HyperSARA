clc; clear all; close all;
format compact;

addpath ../../../../CubeHelix;
addpath ../../../../export_fig_master;
addpath ../../../mnras_faceted_corrected/final_real_data;
addpath ../../../lib/print/;

%% Only load useful channels, create average cube progressively
% addpath ../../mnras_faceted_corrected/final_real_data
% info        = fitsinfo('xsol_FacetedHyperSARA_4-8GHz_NEW.fits');
% rowend      = info.PrimaryData.Size(1);
% colend      = info.PrimaryData.Size(2);
% nchannels   = info.PrimaryData.Size(3);
%
% % xclean = zeros(rowend, colend, 3);
% xhs = zeros(rowend, colend, 3);
% x_mean = zeros(rowend,colend);
% for n = 1:nchannels
%     x = fitsread('xsol_FacetedHyperSARA_4-8GHz_NEW.fits','primary',...
%     'Info', info,...
%     'PixelRegion',{[1 1 rowend], [1 1 colend], n});
%     if n == 1
%         xhs(:,:,1) = x;
%     elseif n == nchannels
%         xhs(:,:,2) = x;
%     end
%     x_mean = x_mean + x;
% end
% xhs(:,:,3) = x_mean/nchannels;
% fitswrite(xhs,'x_fhs_reduced.fits');

%%
% info        = fitsinfo('xsol_SARA_4-8GHz.fits');
% rowend      = info.PrimaryData.Size(1);
% colend      = info.PrimaryData.Size(2);
% nchannels   = info.PrimaryData.Size(3);
% xl1 = zeros(rowend, colend, 3);
% x_mean = zeros(rowend, colend);
% for n = 1:nchannels
%     x = fitsread('xsol_SARA_4-8GHz.fits','primary',...
%     'Info', info,...
%     'PixelRegion',{[1 1 rowend], [1 1 colend], 1});
%     if n == 1
%         xl1(:,:,1) = x;
%     elseif n == nchannels
%         xl1(:,:,2) = x;
%     end
%     x_mean = x_mean + x;
% end
% xl1(:,:,3) = x_mean/nchannels;
% fitswrite(xhs,'x_l1_reduced.fits');
% % would need to do the same with the other cubes

%% Parameters
load_images = 1;
load_residuals = 1;
fig_size = [1000, 1000]; % only change wrt Abdullah's code
shift_colorbar = [0, eps, 0, 0]; % + [left, bottom, width, height] to place it where you want
extension = '.pdf';
map_img = cubehelix(2048);

%% Load images
if load_images
    xhs = fitsread('x_fhs_reduced.fits');
    xl1 = fitsread('xl1_reduced.fits');
    xclean = fitsread('xclean_reduced.fits');
end
[N1, N2, c] = size(xhs);

xclean(:, :, 1) = xclean(:, :, 1) * 42.1827;
xclean(:, :, 2) = xclean(:, :, 2) * 8.3285;
xclean(xclean < 0) = 0;
% xclean = flipud(xclean);
% xhs = flipud(xhs);
% xl1 = flipud(xl1);

%% Load residuals
if load_residuals
    rhs = fitsread('r_fhs_reduced.fits');
    rl1 = fitsread('rl1_reduced.fits');
    rclean = fitsread('rclean_reduced.fits');
end

%% Take only the effective window of the image
a1 = 200; % up
a2 = 250; % down
b1 = 220; % left
b2 = 100;

xhs = xhs(a1:end - a2, b1:end - b2, :);
xl1 = xl1(a1:end - a2, b1:end - b2, :);
xclean = xclean(a1:end - a2, b1:end - b2, :);

rhs = rhs(a1:end - a2, b1:end - b2, :);
rl1 = rl1(a1:end - a2, b1:end - b2, :);
rclean = rclean(a1:end - a2, b1:end - b2, :);

%% Diplay results (first band, last band, average image)

% =========================================================================%
% Plot parameters
% =========================================================================%
mkdir figs;
psf_flux = [40.7 7.97 1];
%% Plot full images
clim_log = [5e-6 0.003   % band 1
            1e-6 0.01    % band end
            1e-6 0.003]; % band average
% clim_log=[5e-6 0.02]; % HS and SARA
% clim_log_cl=psf_flux.*[5e-6 0.02]; % CLEAN
fontsize = 20;

% faceted HyperSARA, SARA
for band = 1:size(xhs, 3)
    [f, h] = display_real_images(xhs(:, :, band), fig_size, shift_colorbar, clim_log(band, :), map_img, fontsize);
    pause(0.5);
    set(h, 'XTick', [1e-5 1e-4 1e-3]);
    export_fig(['figs/xhs_ch', num2str(band), extension], '-transparent', '-q101');
    close;

    [f, h] = display_real_images(xl1(:, :, band), fig_size, shift_colorbar, clim_log(band, :), map_img, fontsize);
    pause(0.5);
    set(h, 'XTick', [1e-5 1e-4 1e-3]);
    export_fig(['figs/xl1_ch', num2str(band), extension], '-transparent', '-q101');
    close;
end

% CLEAN
clim_log = [5e-6 0.02   % band 1
            5e-6 0.02   % band end
            1e-6 1e-2]; % band average
% clim_log=[5e-6 0.02]; % HS and SARA
% clim_log_cl=psf_flux.*[5e-6 0.02]; % CLEAN

for band = 1:size(xhs, 3)
    [f, h] = display_real_images(abs(xclean(:, :, band)), fig_size, shift_colorbar, psf_flux(band) .* clim_log(band, :), map_img, fontsize);
    pause(0.5);
    if band == 3
        set(h, 'XTick', [1e-5 1e-4 1e-3]);
    else
        set(h, 'XTick', [1e-4 1e-3 1e-2 1e-1]);
    end
    export_fig(['figs/xclean_ch', num2str(band), extension], '-transparent', '-q101');
    close;
end

%% Plot residuals

clim_log = [-0.011, 0.011     % band 1
            -0.0055, 0.0055   % band end
            -0.0011, 0.0011]; % band average
fontsize = 40;

for band = 1:size(xhs, 3)
    % faceted hypersara
    [f, h] = display_real_residual(rhs(:, :, band), fig_size, clim_log(band, :), map_img, fontsize);
    if band == 1
        set(h, 'XTick', [-0.01 0 0.01]);
        h.Ruler.Exponent = -2;
    elseif band == 2
        set(h, ...
            'XTick', [-5e-3 0 5e-3]);
    elseif band == 3
        set(h, ...
            'XTick', [-1e-3 0 1e-3]);
    end
    export_fig(['figs/rhs_ch', num2str(band), extension], '-transparent', '-q101');
    close;

    % sara
    [f, h] = display_real_residual(rl1(:, :, band), fig_size, clim_log(band, :), map_img, fontsize);
    if band == 1
        set(h, 'XTick', [-0.01 0 0.01]);
        h.Ruler.Exponent = -2;
    elseif band == 2
        set(h, ...
            'XTick', [-5e-3 0 5e-3]);
    elseif band == 3
        set(h, ...
            'XTick', [-1e-3 0 1e-3]);
    end
    export_fig(['figs/rl1_ch', num2str(band), extension], '-transparent', '-q101');
    close;

    % clean
    [f, h] = display_real_residual(rclean(:, :, band), fig_size, clim_log(band, :), map_img, fontsize);
    if band == 1
        set(h, 'XTick', [-0.01 0 0.01]);
        h.Ruler.Exponent = -2;
    elseif band == 2
        set(h, ...
            'XTick', [-5e-3 0 5e-3]);
    elseif band == 3
        set(h, ...
            'XTick', [-1e-3 0 1e-3]);
    end
    export_fig(['figs/rclean_ch', num2str(band), extension], '-transparent', '-q101');
    close;
end
