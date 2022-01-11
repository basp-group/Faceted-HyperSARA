clc; clear; close all;
format compact;

addpath  '/Users/ad33/CodesScience/Faceted-Hyper-SARA/final_results/CubeHelix';
addpath /Users/ad33/CodesScience/Faceted-Hyper-SARA/final_results;
addpath /Users/ad33/CodesScience/Faceted-Hyper-SARA/final_results/exp1/fhs;
addpath /Users/ad33/CodesScience/Faceted-Hyper-SARA/final_results/exp1/sara;

addpath  /Users/ad33/CodesScience/Faceted-Hyper-SARA/lib/print;

%% Parameters
load_images = 1;
load_residuals = 1;
fig_size = [1000, 1000]; % only change wrt Abdullah's code
shift_colorbar = [0, -1e-2, 0, 0]; % + [left, bottom, width, height] to place it where you want
extension = '.pdf';
extension_png = '.png';

map_img = cubehelix(2048);

%% Take only the effective window of the image
xfhs = [];  xsara = [];
xfhso = fitsread('full_fhs_cube_480.fits');
[N1, N2, c] = size(xfhso);
idFreq = [241 480];
cr = [N1 * 0.5 + 1 N2 * 0.5 + 1];
pixelSZ = 0.06;
fov = [0.06 * 1200; 0.06 * 2500];
fovpix = floor(fov / pixelSZ / 2);
xfhs(:, :, 1) = xfhso(cr(1) - fovpix(1):cr(1) + fovpix(1), cr(2) - fovpix(2):cr(2) + fovpix(2), idFreq(1));
xfhs(:, :, 2) = xfhso(cr(1) - fovpix(1):cr(1) + fovpix(1), cr(2) - fovpix(2):cr(2) + fovpix(2), idFreq(2));
clear xfhso;
xsarao = fitsread('full_sara_cube_480.fits');

xsara(:, :, 1) = xsarao(cr(1) - fovpix(1):cr(1) + fovpix(1), cr(2) - fovpix(2):cr(2) + fovpix(2), idFreq(1));
xsara(:, :, 2) = xsarao(cr(1) - fovpix(1):cr(1) + fovpix(1), cr(2) - fovpix(2):cr(2) + fovpix(2), idFreq(2));
clear xsarao;
% Diplay results (first band, last band, average image)

%% =========================================================================%
% Plot parameters
% =========================================================================%
mkdir figs;
% Plot full images
clim_log = [1e-6 0.002   % band 1
    1e-6 0.002];  % band end
fontsize = 20;
% faceted HyperSARA, SARA
for band = 1:size(xfhs, 3)
    [f, h] = display_image(flipud(xfhs(:, :, band)), fig_size, shift_colorbar, clim_log(band, :), map_img, fontsize, true, 'southoutside');
    set(h, 'XTick', [1e-5 1e-4 1e-3]);
    exportgraphics(f, ['figs/exp1xfhs_ch', num2str(idFreq(band)), extension], 'ContentType', 'vector', ...
        'BackgroundColor', 'none');
    close;
end

for band = 1:2
    [f, h] = display_image(flipud(abs(xsara(:, :, band))), fig_size, shift_colorbar, clim_log(band, :), map_img, fontsize, true, 'southoutside');
    set(h, 'XTick', [1e-5 1e-4 1e-3]);
    exportgraphics(f, ['figs/exp1xsara_ch', num2str(idFreq(band)), extension], 'ContentType', 'vector', ...
        'BackgroundColor', 'none');
    close;
end

%% Plot residuals
%% Plot residuals
clear  xfhso xcleano;

dummy = fitsread('residual_norm_fhs_cube1.fits');
xfhs(:, :, 1) = dummy(cr(1) - fovpix(1):cr(1) + fovpix(1), cr(2) - fovpix(2):cr(2) + fovpix(2), 16);
dummy = fitsread('residual_norm_fhs_cube16.fits');
xfhs(:, :, 2) = dummy(cr(1) - fovpix(1):cr(1) + fovpix(1), cr(2) - fovpix(2):cr(2) + fovpix(2), end);

dummy = fitsread('residual_norm_sara_cube1_ch16.fits');
xsara(:, :, 1) = dummy(cr(1) - fovpix(1):cr(1) + fovpix(1), cr(2) - fovpix(2):cr(2) + fovpix(2));
dummy = fitsread('residual_norm_sara_cube16_ch32.fits');
xsara(:, :, 2) = dummy(cr(1) - fovpix(1):cr(1) + fovpix(1), cr(2) - fovpix(2):cr(2) + fovpix(2));

clim_lin = [-0.005 0.005   % band 1
            -0.002 0.002];  % band end
           % 1e-6 0.003]; % band average
% clim_log=[5e-6 0.02]; % HS and SARA
% clim_log_cl=psf_flux.*[5e-6 0.02]; % CLEAN
fontsize = 20;

% faceted HyperSARA, SARA
for band = 1:size(xfhs, 3)
    'fhs';
    dummy = xfhs(:, :, band);
    std(dummy(:));
    [f, h] = display_image(flipud(xfhs(:, :, band)), fig_size, shift_colorbar, clim_lin(band, :), map_img, fontsize, false, 'southoutside');
    set(h, 'XTick',  [-5e-3 -1e-3 0 1e-3 5e-3]);
    exportgraphics(f, ['figs/res_xhs_ch', num2str(band), extension], 'ContentType', 'vector', ...
    'BackgroundColor', 'none');
  exportgraphics(f, ['figs/res_xhs_ch', num2str(band), extension], 'ContentType', 'vector', ...
    'BackgroundColor', 'none');
    close;
end

% faceted HyperSARA, SARA
for band = 1:size(xfhs, 3)
    'sara';
    dummy = xsara(:, :, band);
    std(dummy(:));
    [f, h] = display_image(flipud(xsara(:, :, band)), fig_size, shift_colorbar, clim_lin(band, :), map_img, fontsize, false, 'southoutside');
    set(h, 'XTick',  [-5e-3 -1e-3 0 1e-3 5e-3]);
    exportgraphics(f, ['figs/res_sara_ch', num2str(band), extension], 'ContentType', 'vector', ...
    'BackgroundColor', 'none');
  exportgraphics(f, ['figs/res_sara_ch', num2str(band), extension], 'ContentType', 'vector', ...
    'BackgroundColor', 'none');
    close;
end

% clim_log = [-0.011, 0.011     % band 1
%             -0.0055, 0.0055   % band end
%
