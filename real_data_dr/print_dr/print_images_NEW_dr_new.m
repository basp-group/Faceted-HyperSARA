clc; clear all; close all;
format compact;

addpath ../../../CubeHelix
addpath ../../../export_fig_master

%% Parameters
load_images = 1;
load_residuals = 1;
fig_size = [1000, 1000]; % only change wrt Abdullah's code
shift_colorbar = [0,-1.5e-2,0,0]; % + [left, bottom, width, height] to place it where you want
extension = '.pdf';

% two regions with the zooms

% inner core
cen0 = [1110 422];
len0 = 120;

% west hotspot
cen = [764 1298];
len = 57;

%% Load images
if load_images
    xhs = fitsread('x_fhs_reduced_dr.fits');
    xl1 = fitsread('xl1_reduced_dr.fits');
    xclean = fitsread('xclean_reduced_dr.fits');
end
[N1, N2, c] = size(xhs);

% xclean(:,:,1) = xclean(:,:,1) * 42.1827;
% xclean(:,:,2) = xclean(:,:,2) * 8.3285;

%% Load residuals
if load_residuals
    rhs = fitsread('r_fhs_reduced_dr.fits');
    rl1 = fitsread('rl1_reduced_dr.fits');
    rclean = fitsread('rclean_reduced_dr.fits');
end

%% Take only the effective window of the image
a1 = 200; %up
a2 = 250; %down
b1 = 220; %left
b2 = 100;

xhs1 = xhs(a1:end-a2,b1:end-b2,:);
xl11 = xl1(a1:end-a2,b1:end-b2,:);
xclean1 = xclean(a1:end-a2,b1:end-b2,:);

rhs1 = rhs(a1:end-a2,b1:end-b2,:);
rl11 = rl1(a1:end-a2,b1:end-b2,:);
rclean1 = rclean(a1:end-a2,b1:end-b2,:);

%% Diplay reasults (images, just first and last band)
% ground truth, faceted hyperSARA, hyperSARA, SARA

%=========================================================================%
% Plot parameters
%=========================================================================%
mkdir figs
% psf_flux = [40.7 7.97 1];
load('flux_psf.mat')
% psf_flux = [40.7 7.97 1];
psf_flux = [flux(1:2), 1];

%% Plot full images
clim_log = [5e-6 0.003;  %band 1
    1e-6 0.01;   %band end
    1e-6 0.003]; %band average

%clim_log=[5e-6 0.02]; % HS and SARA
% clim_log_cl=psf_flux.*[5e-6 0.02]; % CLEAN

map_img = cubehelix(2048);
fontsize=20;

for band = 1:size(xhs,3)
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    imagesc((xhs1(:,:,band)), clim_log(band,:));
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    set(gca,'ColorScale','log')
    h = colorbar;
    set(h,'Location','southoutside');
    set(h,'Fontsize',fontsize);
    if band ==3
        set(h,'XTick',[1e-5 1e-4 1e-3]);
    else
        set(h,'XTick',[1e-5 1e-4 1e-3]);
    end
    pause(0.5)
    h.Position = h.Position + shift_colorbar;
%     aspect = get(ax,'PlotBoxAspectRatio');
%     set(ax,'Units','pixels');
%     pos = get(ax,'Position');
%     pos(3) = aspect(1)/aspect(2)*pos(4);
%     set(ax,'Position',pos);
    % to be corrected!!
%     rectangle('Position',[cen-floor(len/2), len, len],'EdgeColor',[1 1 1],...
%     'LineWidth',2)
%     rectangle('Position',[cen0-floor(len0/2), len0, len0],'EdgeColor',[1 1 1],...
%     'LineWidth',2)
%     rectangle('Position',[0,0, len0, len0],'EdgeColor',[1 1 1],...
%     'LineWidth',2)
    export_fig(['figs/xhs_ch', num2str(band),extension], '-transparent', '-q101')
    
    %
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    imagesc((xl11(:,:,band)), clim_log(band,:));
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Location','southoutside');
    set(h,'Fontsize',fontsize)
    set(gca,'ColorScale','log')
    if band ==3
        set(h,'XTick',[1e-5 1e-4 1e-3]);
    else
        set(h,'XTick',[1e-5 1e-4 1e-3]);
    end
    pause(0.5)
    h.Position = h.Position + shift_colorbar;
    export_fig(['figs/xl1_ch', num2str(band),extension], '-transparent', '-q101')
end
close all

%% Plot full images CLEAN
%% Plot full images
% clim_log = [5e-6 0.02;  %band 1
%             5e-6 0.02;  %band end
%             1e-6 1e-2]; %band average
clim_log = [1e-6 0.0003;  %band 1
    1e-6 0.0003;  %band end
    1e-6 0.003]; %band average

%clim_log=[5e-6 0.02]; % HS and SARA
% clim_log_cl=psf_flux.*[5e-6 0.02]; % CLEAN

map_img = cubehelix(2048);
fontsize=20;

for band = 1:size(xhs,3)
    x1 = xclean1(:,:,band);
    x1(x1 < 0) = 0;
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    imagesc((abs(x1)), psf_flux(band) .* clim_log(band,:));
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Location','southoutside');
    set(h,'Fontsize',fontsize)
    set(gca,'ColorScale','log')
    if band ==3
        set(h,'XTick',[1e-5 1e-4 1e-3]);
    else
        set(h,'XTick',[1e-4 1e-3 1e-2 1e-1]);
    end
    pause(0.5)
    h.Position = h.Position + shift_colorbar;
    export_fig(['figs/xclean_ch', num2str(band),extension], '-transparent', '-q101')
end
close all


%%
%% Plot residuals

clim_log = [-0.011,0.011; %band 1
    -0.0055,0.0055;  %band end
    -0.0011,0.0011]; %band average

map_img = cubehelix(2048);
fontsize=40;

for band = 1:size(xhs,3)
    rhs_z = rhs(:,:,band);
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    im=imagesc((rhs_z), clim_log(band,:));
    colormap(gca, map_img);
    axis image
    ax = gca;
    b = get(ax);
    b =  b.Position; %gets the positon and size of the color bar
    set(gca,'Position',b+[-0.09 -0.026 0.01 0.055])% To change size
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    set(h,'color','white')
    if band ==1
        set(h,'XTick',[-1e-2 0 1e-2])
        h.Ruler.Exponent = -2;
    elseif band==2
        set(h,...
            'XTick',[-5e-3 0 5e-3])
    elseif band==3
        set(h,...
            'XTick',[-1e-3 0 1e-3])
    end
    % Create bounding box around imagesc
    pause(0.5)
    aspect = get(ax,'PlotBoxAspectRatio');
    set(ax,'Units','pixels');
    pos = get(ax,'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(ax,'Position',pos);
    set(ax,'Units','normalized');
    pos = get(ax,'Position');
    annotation(f,'rectangle',...
    [pos(1), pos(2), im.Parent.InnerPosition(3), pos(4)],...
    'Color',[1 1 1],...
    'LineWidth',4);
    export_fig(['figs/rhs_ch', num2str(band),extension],'-transparent', '-q101')
    
    %
    rl1_z = rl1(:,:,band);
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    im=imagesc((rl1_z), clim_log(band,:));
    colormap(gca, map_img);
    axis image
    ax = gca;
    b = get(ax);
    b =  b.Position; %gets the positon and size of the color bar
    set(gca,'Position',b+[-0.09 -0.026 0.01 0.055])% To change size
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    set(h,'color','white')
    if band ==1
        set(h,'XTick',[-1e-2 0 1e-2])
        h.Ruler.Exponent = -2;
    elseif band==2
        set(h,...
            'XTick',[-5e-3 0 5e-3])
    elseif band==3
        set(h,...
            'XTick',[-1e-3 0 1e-3])
    end
    % Create bounding box around imagesc
    pause(0.5)
    aspect = get(ax,'PlotBoxAspectRatio');
    set(ax,'Units','pixels');
    pos = get(ax,'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(ax,'Position',pos);
    set(ax,'Units','normalized');
    pos = get(ax,'Position');
    annotation(f,'rectangle',...
    [pos(1), pos(2), im.Parent.InnerPosition(3), pos(4)],...
    'Color',[1 1 1],...
    'LineWidth',4);
    export_fig(['figs/rl1_ch', num2str(band),extension], '-transparent', '-q101')
    
    %
    rclean_z = rclean(:,:,band);
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    im=imagesc((rclean_z), clim_log(band,:));
    colormap(gca, map_img);
    axis image
    ax = gca;
    b = get(ax);
    b =  b.Position; %gets the positon and size of the color bar
    set(gca,'Position',b+[-0.09 -0.026 0.01 0.055])% To change size
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    set(h,'color','white')
    if band ==1
        set(h,'XTick',[-1e-2 0 1e-2])
        h.Ruler.Exponent = -2;
    elseif band==2
        set(h,...
            'XTick',[-5e-3 0 5e-3])
    elseif band==3
        set(h,...
            'XTick',[-1e-3 0 1e-3])
    end
    % Create bounding box around imagesc
    pause(0.5)
    aspect = get(ax,'PlotBoxAspectRatio');
    set(ax,'Units','pixels');
    pos = get(ax,'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(ax,'Position',pos);
    set(ax,'Units','normalized');
    pos = get(ax,'Position');
    annotation(f,'rectangle',...
    [pos(1), pos(2), im.Parent.InnerPosition(3), pos(4)],...
    'Color',[1 1 1],...
    'LineWidth',4);
    export_fig(['figs/rclean_ch', num2str(band),extension], '-transparent', '-q101')
    
    %waitforbuttonpress
end
close all
