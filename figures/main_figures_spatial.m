clc; clear all; close all;
format compact;

addpath ../../CubeHelix
addpath ../data
addpath ../lib/print/
mkdir figs_spatial
addpath ../results/final/spatial_new
addpath ../lib/faceted-wavelet-transform/src

%% Define experiment type
Ny = 1024;
Nx = 2048;
nChannels = 20;
Qc = 1;
ind = 1;
exp_type = 'spatial';
image_name = 'cygASband_Cube_1024_2048_20';
location = 'southoutside'; % location of the colorbar in the image

%% Load upsilon value
load(strcat('Anorm_hs', ...
    '_Ny=',num2str(Ny),'_Nx=',num2str(Nx), '_L=',num2str(nChannels), ...
    '.mat'), 'operator_norm')
load(strcat('y_', ...
    exp_type,'_',image_name, '_srf=2', ...
    '_Ny=',num2str(Ny),'_Nx=',num2str(Nx),'_L=', num2str(nChannels), ...
    '_snr=40.mat'), 'sigma_noise');
upsilon = (sigma_noise([1,end]).^2)./operator_norm([1,end]); % ! operator_norm already squared

%% Define tile limits
Qx = 4;
Qy = 4;
rg_x = split_range(Qx, Nx);
rg_y = split_range(Qy, Ny);
mark_length = 30; % mark length of the mark
mark_width = 2; % mark length of the mark

%% Parameters
load_images = 1;
load_residuals = 1;
fig_size = [600, 600]; % only change wrt Abdullah's code
shift_colorbar = [0,-1e-2,0,0]; % + [left, bottom, width, height] to place it where you want
extension = '.pdf';
map_img = cubehelix(2048);

img_filenames = {[image_name,'.fits'], 'x_sara.fits', 'x_fhs_Q=4_ovl=0.1.fits', 'x_hs.fits', 'x_fhs_Q=4_ovl=0.fits'};
name = {'x', 'x_l1', 'x_fhs', 'x_hs', 'x_fhs0'};
clim_log = [1e-7 0.5; % band 1
            1e-7 0.5]; % last band 
% clim_log = [1e-5 0.5; % band 1
%             1e-5 0.5]; % last band 
% clim_log = [1e-5 1.2; % band 1
%             1e-5 1.2]; % last band 
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
    max_x = max(x(:));
    
    % insert NaN to create the mark at the borders of the panel
    % ! NaN appear in black in cubehelix !
    for q = 1:Qy-1
        x(rg_y(q, 2)-mark_width:rg_y(q, 2)+mark_width, 1:mark_length, :) = max_x;
        x(rg_y(q, 2)-mark_width:rg_y(q, 2)+mark_width, end-mark_length+1:end, :) = max_x;
    end
    
    for q = 1:Qx-1
        x(1:mark_length, rg_x(q, 2)-mark_width:rg_x(q, 2)+mark_width, :) = max_x;
        x(end-mark_length+1:end, rg_x(q, 2)-mark_width:rg_x(q, 2)+mark_width, :) = max_x;
    end

    % ! see if this is fine
%     if ~strcmp(img_filenames{k}, [image_name,'.fits'])
%         % use ds9 display, using 1/upsilon as a
%         % y = log10(ax+1)/log10(a)
%         x = log(1 + x./reshape(upsilon, [1, 1, 2]))./log(reshape(1./upsilon, [1, 1, 2])); % need to deactivate log colormap when using this
%     end

    for band = 1:size(x, 3)        
        [f, h] = display_image(x(:,:,band), fig_size, shift_colorbar, clim_log(band,:), map_img, fontsize, true, location);
        
        % https://www.mathworks.com/help/matlab/creating_plots/save-figure-at-specific-size-and-resolution.html
        exportgraphics(f,['figs_spatial/', name{k}, num2str(band), extension],'ContentType','vector',...
               'BackgroundColor','none')
        close
    end
end
         
%% Plot residuals
res_filenames = {'res_sara.fits', 'res_fhs_Q=4_ovl=0.1.fits', 'res_hs.fits', 'res_fhs_Q=4_ovl=0.fits'};
name = {'res_l1', 'res_fhs', 'res_hs', 'res_fhs0'};
clim_log = [-1e-5,1e-5;     %band 1
            -1e-5,1e-5];  %band end
% fig_size = [1000, 1000];
% fontsize = 40;
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
        [f, h] = display_image(res(:,:,band), fig_size, shift_colorbar, clim_log(band,:), map_img, fontsize, false, location);
        % if band ==1
        %     set(h,'XTick',[-0.01 0 0.01])
        %     h.Ruler.Exponent = -2;
        % elseif band==2
        %     set(h,...
        %         'XTick',[-5e-3 0 5e-3])
        % elseif band==3
        %     set(h,...
        %         'XTick',[-1e-3 0 1e-3])
        % end
        
        % https://www.mathworks.com/help/matlab/creating_plots/save-figure-at-specific-size-and-resolution.html
        exportgraphics(f,['figs_spatial/', name{k}, num2str(band),extension],'ContentType','vector',...
               'BackgroundColor','none')
        close

    end
end
