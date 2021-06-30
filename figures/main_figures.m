clc; clear all; close all;
format compact;

addpath ../../CubeHelix
addpath ../../export_fig_master
addpath ../data
addpath ../lib/print/
mkdir figs
addpath ../results/spatial

%% Define experiment type
Ny = 1024;
Nx = 2048;
nChannels = 20;
Qc = 1;
ind = 1;
exp_type = 'spatial';
image_name = 'cygASband_Cube_1024_2048_20';

%% Load upsilon value
load(strcat('Anorm_hs', ...
    '_Ny=',num2str(Ny),'_Nx=',num2str(Nx), '_L=',num2str(nChannels), ...
    '_Qc=',num2str(Qc),'_ind=',num2str(ind), '_ch=', num2str(ind), '.mat'), 'operator_norm')
load(strcat('y_', ...
    exp_type,'_',image_name, '_srf=2', ...
    '_Ny=',num2str(Ny),'_Nx=',num2str(Nx),'_L=', num2str(nchannels), ...
    '_snr=40.mat'), 'sigma_noise');
upsilon = (sigma_noise([1,end]).^2)./operator_norm([1,end]); % ! operator_norm already squared

%% Parameters
load_images = 1;
load_residuals = 1;
fig_size = [1000, 1000]; % only change wrt Abdullah's code
shift_colorbar = [0,eps,0,0]; % + [left, bottom, width, height] to place it where you want
extension = '.pdf';
map_img = cubehelix(2048);

img_filenames = {[image_name,'.fits'], 'x_sara.fits', 'x_fhs.fits', 'x_hs.fits'};
name = {'x', 'x_l1', 'x_fhs', 'x_hs'};
clim_log = [5e-6 0.003; % band 1
            1e-6 0.01]; % last band 
fontsize=20;
%% 

for k = 1:numel(img_filename)
    info        = fitsinfo(img_filenames{k});
    rowend      = info.PrimaryData.Size(1);
    colend      = info.PrimaryData.Size(2);
    nchannels   = info.PrimaryData.Size(3);

    x = fitsread(img_filenames{k},'primary',...
        'Info', info,...
        'PixelRegion',{[1 1 rowend], [1 1 colend], [1, colend-1, colend]});

    % ! see if this is fine
    x = log10(1 + x./reshape(upsilon, [1, 1, 2]));    

    for band = 1:size(x, 3)        
        [f, h] = display_image2(x(:,:,band), fig_size, shift_colorbar, clim_log(band,:), map_img, fontsize);
        % pause(0.5)
        % set(h,'XTick',[1e-5 1e-4 1e-3]);
        export_fig(['figs/', name{k}, num2str(band), extension], '-transparent', '-q101')
        close
    end
end
         
%% Plot residuals
res_filenames = {'res_sara.fits', 'res_fhs.fits', 'res_hs.fits'};
name = {'res_l1', 'res_fhs', 'res_hs'};
clim_log = [-0.011,0.011;    %band 1
            -0.0055,0.0055];  %band end
fontsize = 40;

for k = 1:numel(res_filename)
    info        = fitsinfo(res_filenames{k});
    rowend      = info.PrimaryData.Size(1);
    colend      = info.PrimaryData.Size(2);
    nchannels   = info.PrimaryData.Size(3);

    res = fitsread(res_filenames{k},'primary',...
        'Info', info,...
        'PixelRegion',{[1 1 rowend], [1 1 colend], [1, colend-1, colend]});

    for band = 1:size(xhs,3)
        % faceted hypersara
        [f, h] = display_real_residual(res(:,:,band), fig_size, clim_log(band,:), map_img, fontsize);
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
        export_fig(['figs/', name{k}, num2str(band),extension],'-transparent', '-q101')
        close

    end
end
