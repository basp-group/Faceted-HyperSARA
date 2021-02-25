clc; clear all; close all;
format compact;

addpath ../../../CubeHelix
addpath ../../../export_fig_master

%% Parameters
load_images = 1;
load_residuals = 1;
fig_size = [1250, 1250]; % only change wrt Abdullah's code
shift_colorbar = [0,0,0,0]; % + [left, bottom, width, height] to place it where you want
pixel_shift_colorbar = 30;
% [0,-2.5e-2,0,0]
% [-0.065 -0.026 0.015 0.055]
extension = '.pdf';

%% Load images
if load_images
    xhs = fitsread('x_fhs_reduced_dr.fits');
    xl1 = fitsread('xl1_reduced_dr.fits');
    xclean = fitsread('xclean_reduced_dr.fits');
end
[N1, N2, c] = size(xhs);

xclean(:,:,1) = xclean(:,:,1) * 42.1827;
xclean(:,:,2) = xclean(:,:,2) * 8.3285;

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
psf_flux = [flux(1:2), 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% skip zooms for now %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the east hot spot
clim_log = [0.0005,0.14]; % HS and SARA
tickLabel = [0.001,0.01,0.1];
%clim_log_cl = psf_flux * [0.008 0.05]; % CLEAN
%tickLabel_cl = psf_flux(band) * [0.01,0.04];

map_img = cubehelix(2048);
fontsize=60;
cen = [410 2282];
len = 130;

for band = 1:size(xhs,3)
    xhs_z = xhs(:,:,band);
    xhs_z1 = xhs_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    im = imagesc((xhs_z1), clim_log);
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    set(h,'color','white')
    set(gca,'ColorScale','log');
    set(h,'XTick',tickLabel);
    % a=get(h);
    % a =  a.Position; %gets the positon and size of the color bar
    % set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
    % Create rectangle
%     annotation(f,'rectangle',...
%     [0.2466 0.106995174129353 0.539200000000001 0.82089552238806],...
%     'Color',[1 1 1],...
%     'LineWidth',4);    
    pause(0.5)
    % Get the new aspect ratio data
    aspect = get(ax,'PlotBoxAspectRatio');
    % Change axes Units property (this only works with non-normalized units)
    set(ax,'Units','pixels');
    % Get and update the axes width to match the plot aspect ratio.
    pos = get(ax,'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(ax,'Position',pos);
    set(h,'Units', 'pixels')
    pos_h = get(h,'Position');
    pos_h = [pos(1)+pos(3)+pixel_shift_colorbar, pos_h(2), 2*pos_h(3), pos(4)];
    set(h,'Position',pos_h);
    set(ax,'Units','normalized');
    pos = get(ax,'Position');
    annotation(f,'rectangle',...
    [pos(1), pos(2), im.Parent.InnerPosition(3), pos(4)],...
    'Color',[1 1 1],...
    'LineWidth',4);
    export_fig(['figs/xhs_z11_ch', num2str(band),extension],'-transparent','-q101')
    
    %
    xl1_z = xl1(:,:,band);
    xl1_z1 = xl1_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    im=imagesc((xl1_z1), clim_log);
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    set(h,'color','white')
    set(gca,'ColorScale','log');
    set(h,'XTick',tickLabel);
    % a=get(h);
    % a =  a.Position; %gets the positon and size of the color bar
    % set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
    % Create bounding box around imagesc
    pause(0.5)
    aspect = get(ax,'PlotBoxAspectRatio');
    set(ax,'Units','pixels');
    pos = get(ax,'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(ax,'Position',pos);
    set(h,'Units', 'pixels')
    pos_h = get(h,'Position');
    pos_h = [pos(1)+pos(3)+pixel_shift_colorbar, pos_h(2), 2*pos_h(3), pos(4)];
    set(h,'Position',pos_h);
    set(ax,'Units','normalized');
    pos = get(ax,'Position');
    annotation(f,'rectangle',...
    [pos(1), pos(2), im.Parent.InnerPosition(3), pos(4)],...
    'Color',[1 1 1],...
    'LineWidth',4);
    export_fig(['figs/xl1_z11_ch', num2str(band),extension], '-transparent','-q101')
    
    %
%     xclean_z = xclean(:,:,band);
%     xclean_z(xclean_z < 0) = 0;
%     xclean_z1 = xclean_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
%     f=figure('visible','on');
%     set(gca, 'Color', 'none'); % sets axes background
%     set(f,'PaperUnits','centimeters')
%     set(f,'PaperType','A4');
%     set(f,'PaperOrientation',orient);
%     set(f,'units','pixel','outerposition',[0 0 fig_size])
%     imagesc((xclean_z1), psf_flux(band) .* clim_log);
%     colormap(gca, map_img);
%     axis image
%     ax = gca;
%     ax.YAxis.Visible = 'off';
%     ax.XAxis.Visible = 'off';
%     ax.TickLabelInterpreter='latex';
%     h = colorbar;
%     set(h,'Fontsize',fontsize)
%     set(h,'color','white')
%     set(gca,'ColorScale','log');
%     set(h,'XTick',round(psf_flux(band) .* [0.01,0.02,0.03,0.04],2));
%     a=get(h);
%     a =  a.Position; %gets the positon and size of the color bar
%     set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
%     % Create rectangle
%     annotation(f,'rectangle',...
%     [0.2466 0.106995174129353 0.539200000000001 0.82089552238806],...
%     'Color',[1 1 1],...
%     'LineWidth',4);
%     export_fig(['figs/xclean_z1_ch', num2str(band),extension], '-transparent','-q101')
    
    %waitforbuttonpress
end
close all

%% Plot the west hot spot
clim_log = [0.0005,0.2]; % HS and SARA
tickLabel = [0.001,0.01,0.1];
%clim_log_cl = psf_flux * [0.008 0.05]; % CLEAN
%tickLabel_cl = psf_flux(band) * [0.01,0.04];

map_img = cubehelix(2048);
fontsize=60;
cen = [1110 422];
len = 120;

for band = 1:size(xhs,3)
    xhs_z = xhs(:,:,band);
    xhs_z1 = xhs_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    im=imagesc((xhs_z1), clim_log);
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    set(h,'color','white')
    set(gca,'ColorScale','log');
    set(h,'XTick',tickLabel);
    % a=get(h);
    % a =  a.Position; %gets the positon and size of the color bar
    % set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
    % Create bounding box around imagesc
    pause(0.5)
    aspect = get(ax,'PlotBoxAspectRatio');
    set(ax,'Units','pixels');
    pos = get(ax,'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(ax,'Position',pos);
    set(h,'Units', 'pixels')
    pos_h = get(h,'Position');
    pos_h = [pos(1)+pos(3)+pixel_shift_colorbar, pos_h(2), 2*pos_h(3), pos(4)];
    set(h,'Position',pos_h);
    set(ax,'Units','normalized');
    pos = get(ax,'Position');
    annotation(f,'rectangle',...
    [pos(1), pos(2), im.Parent.InnerPosition(3), pos(4)],...
    'Color',[1 1 1],...
    'LineWidth',4);
    export_fig(['figs/xhs_z31_ch', num2str(band),extension],'-transparent','-q101')
    
    %
    xl1_z = xl1(:,:,band);
    xl1_z1 = xl1_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    im=imagesc((xl1_z1), clim_log);
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    set(h,'color','white')
    set(gca,'ColorScale','log');
    set(h,'XTick',tickLabel);
    % a=get(h);
    % a =  a.Position; %gets the positon and size of the color bar
    % set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
    % Create bounding box around imagesc
    pause(0.5)
    aspect = get(ax,'PlotBoxAspectRatio');
    set(ax,'Units','pixels');
    pos = get(ax,'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(ax,'Position',pos);
    set(h,'Units', 'pixels')
    pos_h = get(h,'Position');
    pos_h = [pos(1)+pos(3)+pixel_shift_colorbar, pos_h(2), 2*pos_h(3), pos(4)];
    set(h,'Position',pos_h);
    set(ax,'Units','normalized');
    pos = get(ax,'Position');
    annotation(f,'rectangle',...
    [pos(1), pos(2), im.Parent.InnerPosition(3), pos(4)],...
    'Color',[1 1 1],...
    'LineWidth',4);
    export_fig(['figs/xl1_z31_ch', num2str(band),extension], '-transparent','-q101')
    
    %
%     xclean_z = xclean(:,:,band);
%     xclean_z(xclean_z < 0) = 0;
%     xclean_z1 = xclean_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
%     f=figure('visible','on');
%     set(gca, 'Color', 'none'); % sets axes background
%     set(f,'PaperUnits','centimeters')
%     set(f,'PaperType','A4');
%     set(f,'PaperOrientation',orient);
%     set(f,'units','pixel','outerposition',[0 0 fig_size])
%     imagesc((xclean_z1), psf_flux(band) .* clim_log);
%     colormap(gca, map_img);
%     axis image
%     ax = gca;
%     ax.YAxis.Visible = 'off';
%     ax.XAxis.Visible = 'off';
%     ax.TickLabelInterpreter='latex';
%     h = colorbar;
%     set(h,'Fontsize',fontsize)
%     set(h,'color','white')
%     set(gca,'ColorScale','log');
%     set(h,'XTick',round(psf_flux(band) .* [0.01,0.02,0.03,0.04],2));
%     a=get(h);
%     a =  a.Position; %gets the positon and size of the color bar
%     set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
%     % Create rectangle
%     annotation(f,'rectangle',...
%     [0.2466 0.106995174129353 0.539200000000001 0.82089552238806],...
%     'Color',[1 1 1],...
%     'LineWidth',4);
%     export_fig(['figs/xclean_z3_ch', num2str(band),extension], '-transparent','-q101')
    
    %waitforbuttonpress
end
close all

%% Plot the inner core
clim_log = [1e-5  0.02];
tickLabel = [1e-4,1e-3,1e-2];
%clim_log_cl = psf_flux * [0.00005  0.0015]; % clean
%tickLabel_cl = round(psf_flux(band) .*[1e-4,1e-3]);

map_img = cubehelix(2048);
fontsize=60;
cen = [764 1298];
len = 57;

for band = 1 : size(xhs,3)
    xhs_z = xhs(:,:,band);
    xhs_z1 = xhs_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    im=imagesc((xhs_z1), clim_log);
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    set(h,'color','white')
    set(h,'XTick',[1e-4,1e-3]);
    set(gca,'ColorScale','log');
    % a=get(h);
    % a =  a.Position; %gets the positon and size of the color bar
    % set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
    % Create ellipse
%     annotation(f,'ellipse',...
%     [0.550090909090909 0.388149939540508 0.0563090909090906 0.0763624982704378],...
%     'Color',[1 1 1],...
%     'LineWidth',5,...
%     'LineStyle','--');
%     % Create rectangle
    % Create bounding box around imagesc
    pause(0.5)
    aspect = get(ax,'PlotBoxAspectRatio');
    set(ax,'Units','pixels');
    pos = get(ax,'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(ax,'Position',pos);
    set(h,'Units', 'pixels')
    pos_h = get(h,'Position');
    pos_h = [pos(1)+pos(3)+pixel_shift_colorbar, pos_h(2), 2*pos_h(3), pos(4)];
    set(h,'Position',pos_h);
    set(ax,'Units','normalized');
    pos = get(ax,'Position');
    annotation(f,'rectangle',...
    [pos(1), pos(2), im.Parent.InnerPosition(3), pos(4)],...
    'Color',[1 1 1],...
    'LineWidth',4);
    export_fig(['figs/xhs_z21_ch', num2str(band),extension], '-transparent','-q101')
    
    %
    xl1_z = xl1(:,:,band);
    xl1_z1 = xl1_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 fig_size])
    im=imagesc((xl1_z1), clim_log);
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    ax.TickLabelInterpreter='latex';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    set(h,'color','white')
    set(h,'XTick',[1e-4,1e-3]);
    set(gca,'ColorScale','log');
    % a=get(h);
    % a =  a.Position; %gets the positon and size of the color bar
    % set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
    % Create ellipse
%     annotation(f,'ellipse',...
%     [0.550090909090909 0.388149939540508 0.0563090909090906 0.0763624982704378],...
%     'Color',[1 1 1],...
%     'LineWidth',5,...
%     'LineStyle','--');
%     % Create rectangle
    % Create bounding box around imagesc
    pause(0.5)
    aspect = get(ax,'PlotBoxAspectRatio');
    set(ax,'Units','pixels');
    pos = get(ax,'Position');
    pos(3) = aspect(1)/aspect(2)*pos(4);
    set(ax,'Position',pos);
    set(h,'Units', 'pixels')
    pos_h = get(h,'Position');
    pos_h = [pos(1)+pos(3)+pixel_shift_colorbar, pos_h(2), 2*pos_h(3), pos(4)];
    set(h,'Position',pos_h);
    set(ax,'Units','normalized');
    pos = get(ax,'Position');
    annotation(f,'rectangle',...
    [pos(1), pos(2), im.Parent.InnerPosition(3), pos(4)],...
    'Color',[1 1 1],...
    'LineWidth',4);
    export_fig(['figs/xl1_z21_ch', num2str(band),extension], '-transparent','-q101')
end
close all