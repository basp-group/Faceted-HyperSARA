clc; clear all; close all;
format compact;

addpath ../../../CubeHelix
addpath ../../../export_fig_master
addpath ../../lib/print

%%
% addpath ../../mnras_faceted_corrected/final_dr_gamma5e-6
% filename = 'facethyperSARA_xsol_it10752_reweight20_gamma5e-06_gamma0_0.001_2b_fouRed2_perc15_adpteps0.fits';
% info        = fitsinfo(filename);
% rowend      = info.PrimaryData.Size(1);
% colend      = info.PrimaryData.Size(2);
% nchannels   = info.PrimaryData.Size(3);
%
% % xclean = zeros(rowend, colend, 3);
% xhs = zeros(rowend, colend, 3);
% x_mean = zeros(rowend,colend);
% for n = 1:nchannels
%     x = fitsread(filename,'primary',...
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
% fitswrite(xhs,'x_fhs_reduced_dr.fits');
%
% filename = 'facethyperSARA_xsol_it10752_reweight20_gamma5e-06_gamma0_0.001_2b_fouRed2_perc15_adpteps0.fits';
% x = fitsread(filename);
% xhs = zeros(size(x,1),size(x,2),3);
% xhs(:,:,1:2) = x(:,:,[1,end]);
% xhs(:,:,3) = mean(x, 3);

%% Parameters
load_images = 1;
load_residuals = 1;
fig_size = [1000, 1000]; % only change wrt Abdullah's code
shift_colorbar = [0,eps,0,0]; % + [left, bottom, width, height] to place it where you want % -1.5e-2
extension = '.pdf';
map_img = cubehelix(2048);

%% Load images
if load_images
    xhs = fitsread('x_fhs_reduced_dr.fits');
    xhs_avg = fitsread('x_fhs_avg_reduced.fits');
    xl1 = fitsread('xl1_reduced_dr.fits');
    xclean = fitsread('xclean_reduced_dr.fits');
end
xhs = flipud(xhs);
xclean(xclean < 0) = 0;
[N1, N2, c] = size(xhs);

%% Load residuals
if load_residuals
    rhs = fitsread('r_fhs_reduced_dr.fits');
    rhs_avg = fitsread('r_fhs_avg_reduced.fits');
    rl1 = fitsread('rl1_reduced_dr.fits');
    rclean = fitsread('rclean_reduced_dr.fits');
end
rhs = flipud(rhs);

%% Take only the effective window of the image
a1 = 200; %up
a2 = 250; %down
b1 = 220; %left
b2 = 100;

xhs = xhs(a1:end-a2,b1:end-b2,:);
xhs_avg = xhs_avg(a1:end-a2,b1:end-b2,:);
xl1 = xl1(a1:end-a2,b1:end-b2,:);
xclean = xclean(a1:end-a2,b1:end-b2,:);

rhs = rhs(a1:end-a2,b1:end-b2,:);
rhs_avg = rhs_avg(a1:end-a2,b1:end-b2,:);
rl1 = rl1(a1:end-a2,b1:end-b2,:);
rclean = rclean(a1:end-a2,b1:end-b2,:);

%% Diplay results (first band, last band, average image)

%=========================================================================%
% Plot parameters
%=========================================================================%
mkdir figs
load('flux_psf.mat')
psf_flux = [flux(1:2), 1];

%% Plot full images

% SARA family
clim_log = [5e-6 0.003;  %band 1
    1e-6 0.01;   %band end
    1e-6 0.003]; %band average
fontsize=20;

for band = 1:size(xhs,3)
    [f, h] = display_real_images(xhs(:,:,band), fig_size, shift_colorbar, clim_log(band,:), map_img, fontsize);
    pause(0.5)
    set(h,'XTick',[1e-5 1e-4 1e-3]);
    export_fig(['figs/xhs_ch', num2str(band),extension], '-transparent', '-q101')
    close
    
    [f, h] = display_real_images(xhs_avg(:,:,band), fig_size, shift_colorbar, clim_log(band,:), map_img, fontsize);
    pause(0.5)
    set(h,'XTick',[1e-5 1e-4 1e-3]);
    export_fig(['figs/xhs_avg_ch', num2str(band),extension], '-transparent', '-q101')
    close
    
    [f, h] = display_real_images(xl1(:,:,band), fig_size, shift_colorbar, clim_log(band,:), map_img, fontsize);
    pause(0.5)
    set(h,'XTick',[1e-5 1e-4 1e-3]);
    export_fig(['figs/xl1_ch', num2str(band),extension], '-transparent', '-q101')
    close
end

% CLEAN
% clim_log = [5e-6 0.02;  %band 1
%             5e-6 0.02;  %band end
%             1e-6 1e-2]; %band average
clim_log = [1e-6 0.0003;  %band 1
    1e-6 0.0003;  %band end
    1e-6 0.003];  %band average

for band = 1:size(xhs,3)
    [f, h] = display_real_images(xclean(:,:,band), fig_size, shift_colorbar, psf_flux(band).*clim_log(band,:), map_img, fontsize);
    pause(0.5)
    if band ==3
        set(h,'XTick',[1e-5 1e-4 1e-3]);
    else
        set(h,'XTick',[1e-4 1e-3 1e-2 1e-1]);
    end
    export_fig(['figs/xclean_ch', num2str(band),extension], '-transparent', '-q101')
    close
end

%% Plot residuals

% clim_log = [-0.011,0.011; %band 1
%     -0.0055,0.0055;       %band end
%     -0.0011,0.0011];      %band average
clim_log = [-1.1e-2,1.1e-2; %band 1 -> could go to 1e-3 for SARA, faceted HyperSARA
    -3.3e-3,3.3e-3;       %band end
    -1.1e-3,1.1e-3];      %band average
fontsize = 40;

for band = 1:size(xhs,3)
    % faceted hypersara
    [f, h] = display_real_residual(rhs(:,:,band), fig_size, clim_log(band,:), map_img, fontsize);
%     if band == 2
%         set(h,'XTick',[-3e-3, 0, 3e-3]);
%         h.Ruler.Exponent = -3;
%     else
%         set(h,'XTick',[round(clim_log(band,1),round(abs(log10(abs(clim_log(band,2)))))) 0 round(clim_log(band,2),round(abs(log10(clim_log(band,2)))))])
%         h.Ruler.Exponent = log10(clim_log(band,2));
%     end
    if band ==1
        set(h,'XTick',[-0.01 0 0.01])
        h.Ruler.Exponent = -2;
    elseif band==2
        set(h,...
            'XTick',[-3e-3 0 3e-3])
    elseif band==3
        set(h,...
            'XTick',[-1e-3 0 1e-3])
    end
    export_fig(['figs/rhs_ch', num2str(band),extension],'-transparent', '-q101')
    close
    
    % average cube faceted hypersara, no DR
    [f, h] = display_real_residual(rhs_avg(:,:,band), fig_size, clim_log(band,:), map_img, fontsize);
    if band ==1
        set(h,'XTick',[-0.01 0 0.01])
        h.Ruler.Exponent = -2;
    elseif band==2
        set(h,...
            'XTick',[-3e-3 0 3e-3])
    elseif band==3
        set(h,...
            'XTick',[-1e-3 0 1e-3])
    end
    export_fig(['figs/rhs_avg_ch', num2str(band),extension],'-transparent', '-q101')
    close
    
    % sara
    [f, h] = display_real_residual(rl1(:,:,band), fig_size, clim_log(band,:), map_img, fontsize);
    if band ==1
        set(h,'XTick',[-0.01 0 0.01])
        h.Ruler.Exponent = -2;
    elseif band==2
        set(h,...
            'XTick',[-3e-3 0 3e-3])
    elseif band==3
        set(h,...
            'XTick',[-1e-3 0 1e-3])
    end
    export_fig(['figs/rl1_ch', num2str(band),extension], '-transparent', '-q101')
    close
    
    % clean
    [f, h] = display_real_residual(rclean(:,:,band), fig_size, clim_log(band,:), map_img, fontsize);
    if band ==1
        set(h,'XTick',[-0.01 0 0.01])
        h.Ruler.Exponent = -2;
    elseif band==2
        set(h,...
            'XTick',[-3e-3 0 3e-3])
    elseif band==3
        set(h,...
            'XTick',[-1e-3 0 1e-3])
    end
    export_fig(['figs/rclean_ch', num2str(band),extension], '-transparent', '-q101')
    close
end
