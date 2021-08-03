clc; clear all; close all;
format compact;

addpath ../../CubeHelix
addpath ../data
addpath ../lib/print/
mkdir figs_spectral
addpath ../results/final/spectral

%% Define experiment type
Ny = 256;
Nx = 512;
nChannels = 100;
Qc = 1;
ind = 1;
exp_type = 'spectral';
image_name = 'cygASband_Cube_256_512_100';

%% Load upsilon value
load(strcat('Anorm_hs', ...
    '_Ny=',num2str(Ny),'_Nx=',num2str(Nx), '_L=',num2str(nChannels), ...
    '.mat'), 'operator_norm')
load(strcat('y_', ...
    exp_type,'_',image_name, '_srf=2', ...
    '_Ny=',num2str(Ny),'_Nx=',num2str(Nx),'_L=', num2str(nChannels), ...
    '_snr=40.mat'), 'sigma_noise');
upsilon = (sigma_noise([1,end]).^2)./operator_norm([1,end]); % ! operator_norm already squared

%% Parameters
load_images = 1;
load_residuals = 1;
fig_size = [600, 600]; % only change wrt Abdullah's code
shift_colorbar = [0,-1e-2,0,0]; % + [left, bottom, width, height] to place it where you want
extension = '.pdf';
map_img = cubehelix(2048);

img_filenames = {[image_name,'.fits'], 'x_sara.fits', 'x_fhs_Qc=1.fits', 'x_fhs_Qc=2.fits', 'x_fhs_Qc=5.fits', 'x_fhs_Qc=10.fits'};
name = {'x', 'x_l1', 'x_fhs_Qc1', 'x_fhs_Qc2', 'x_fhs_Qc5', 'x_fhs_Qc10'};
clim_log = [1e-5 0.5;  % band 1
            1e-5 0.5]; % last band 
fontsize=20;

%% 

for k = 1:numel(img_filenames)
    info        = fitsinfo(img_filenames{k});
    rowend      = info.PrimaryData.Size(1);
    colend      = info.PrimaryData.Size(2);
    nchannels   = info.PrimaryData.Size(3);

    x = fitsread(img_filenames{k},'primary',...
        'Info', info,...
        'PixelRegion',{[1 1 rowend], [1 1 colend], [1, nchannels-1, nchannels]});
    x = flipud(x);

    % ! see if this is fine
%     if ~strcmp(img_filenames{k}, [image_name,'.fits'])
%         % use ds9 display, using 1/upsilon as a
%         % y = log10(ax+1)/log10(a)
%         x = log(1 + x./reshape(upsilon, [1, 1, 2]))./log(reshape(1./upsilon, [1, 1, 2])); % need to deactivate log colormap when using this
%     end

    for band = 1:size(x, 3)        
        [f, h] = display_image(x(:,:,band), fig_size, shift_colorbar, clim_log(band,:), map_img, fontsize, true);
        
        % https://www.mathworks.com/help/matlab/creating_plots/save-figure-at-specific-size-and-resolution.html
        exportgraphics(f,['figs_spectral/', name{k}, num2str(band), extension],'ContentType','vector',...
               'BackgroundColor','none')
        close
    end
end
         
%% Plot residuals
res_filenames = {'res_sara.fits', 'res_fhs_Qc=1.fits', 'res_fhs_Qc=2.fits', 'res_fhs_Qc=5.fits', 'res_fhs_Qc=10.fits'};
name = {'res_l1', 'res_fhs_Qc1', 'res_fhs_Qc2', 'res_fhs_Qc5', 'res_fhs_Qc10'};
clim_log = [-1e-5,1e-5;   % band 1
            -1e-5,1e-5];  % band end
fig_size = [600, 600];
fontsize = 20;

for k = 1:numel(res_filenames)
    info        = fitsinfo(res_filenames{k});
    rowend      = info.PrimaryData.Size(1);
    colend      = info.PrimaryData.Size(2);
    nchannels   = info.PrimaryData.Size(3);

    res = fitsread(res_filenames{k},'primary',...
        'Info', info,...
        'PixelRegion',{[1 1 rowend], [1 1 colend], [1, nchannels-1, nchannels]});
    res = flipud(res);

    for band = 1:size(res,3)
        % faceted hypersara
        [f, h] = display_image(res(:,:,band), fig_size, shift_colorbar, clim_log(band,:), map_img, fontsize, false);
        
        % https://www.mathworks.com/help/matlab/creating_plots/save-figure-at-specific-size-and-resolution.html
        exportgraphics(f,['figs_spectral/', name{k}, num2str(band),extension],'ContentType','vector',...
               'BackgroundColor','none')
        close

    end
end
