clc; clear all; close all;
format compact;

addpath CubeHelix
addpath export_fig_master

%% Parameters
load_images = 1;
load_residuals = 1;

%% Load images
if load_images
    xhs = fitsread('x_fhs_reduced.fits');
    xl1 = fitsread('xl1_reduced.fits');
    xclean = fitsread('xclean_reduced.fits');
end
[N1, N2, c] = size(xhs);

xclean(:,:,1) = xclean(:,:,1) * 42.1827;
xclean(:,:,2) = xclean(:,:,2) * 8.3285;

%% Load residuals
if load_residuals
    rhs = fitsread('r_fhs_reduced.fits');
    rl1 = fitsread('rl1_reduced.fits');
    rclean = fitsread('rclean_reduced.fits');
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
psf_flux = [40.7 7.97 1];
%% Plot full images
clim_log = [5e-6 0.003;  %band 1
            1e-6 0.01;    %band end
            1e-6 0.003]; %band average
        
%clim_log=[5e-6 0.02]; % HS and SARA
% clim_log_cl=psf_flux.*[5e-6 0.02]; % CLEAN

map_img = cubehelix(2048);
fontsize=20;

for band = 3 %:size(xhs,3)
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 1000 1000])
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
    a=get(h);
    a1 =  a.Position; %gets the positon and size of the color bar
    set(h,'Position',[0.13    0.27   0.7750    0.0164]);% To change size
%     % Create rectangle
%     annotation(f,'rectangle',...
%         [0.181 0.34 0.035 0.0426555252756772],'Color',[1 1 1],...
%         'LineWidth',1);
%     
%     % Create rectangle
%     annotation(f,'rectangle',...
%         [0.494 0.5 0.013 0.0169286577992742],'Color',[1 1 1],...
%         'LineWidth',1);
     export_fig(['figs/xhs_ch', num2str(band),'.png'], '-transparent')
    
    %
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 1000 1000])
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
    %     a=get(h);
    %     a =  a.Position; %gets the positon and size of the color bar
    set(h,'Position',[0.13    0.27   0.7750    0.0164]);
%         % Create rectangle
%     annotation(f,'rectangle',...
%         [0.181 0.34 0.035 0.0426555252756772],'Color',[1 1 1],...
%         'LineWidth',1);
%     
%     % Create rectangle
%     annotation(f,'rectangle',...
%         [0.494 0.5 0.013 0.0169286577992742],'Color',[1 1 1],...
%         'LineWidth',1);
    export_fig(['figs/xl1_ch', num2str(band),'.png'], '-transparent')
    
    
    %waitforbuttonpress;
end

%% Plot full images CLEAN
%% Plot full images
clim_log = [5e-6 0.02;  %band 1
            5e-6 0.02;    %band end
            1e-6 1e-2]; %band average
        
%clim_log=[5e-6 0.02]; % HS and SARA
% clim_log_cl=psf_flux.*[5e-6 0.02]; % CLEAN

map_img = cubehelix(2048);
fontsize=20;

for band = 1:size(xhs,3)
    
    %
    x1 = xclean1(:,:,band);
    x1(x1 < 0) = 0;
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 1000 1000])
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
    a=get(h);
    a =  a.Position; %gets the positon and size of the color bar
    set(h,'Position',[0.13    0.27   0.7750    0.0164]);
%     % Create rectangle
%     annotation(f,'rectangle',...
%         [0.181 0.34 0.035 0.0426555252756772],'Color',[1 1 1],...
%         'LineWidth',1);
%     
%     % Create rectangle
%     annotation(f,'rectangle',...
%         [0.494 0.5 0.013 0.0169286577992742],'Color',[1 1 1],...
%         'LineWidth',1);
     export_fig(['figs/xclean_ch', num2str(band),'.png'], '-transparent')
    
    %waitforbuttonpress;
end
         

%%
%% Plot residuals

clim_log = [-0.011,0.011;  %band 1
            -0.0055,0.0055;    %band end
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
    set(f,'units','pixel','outerposition',[0 0 1000 1000])
    imagesc((rhs_z), clim_log(band,:));
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
    a=get(h);
    a =  a.Position; %gets the positon and size of the color bar
    set(h,'Position',a+[0.05 -0.055 0.008 0.113])% To change size
    if band ==1
        set(h,'XTick',[-0.01 0 0.01])
    elseif band==2
        set(h,...
            'XTick',[-5e-3 0 5e-3])
    elseif band==3
        set(h,...
            'XTick',[-1e-3 0 1e-3])
    end
    export_fig(['figs/rhs_ch', num2str(band),'.png'],'-transparent')
    
    %
    rl1_z = rl1(:,:,band);
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 1000 1000])
    imagesc((rl1_z), clim_log(band,:));
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
    a=get(h);
    a =  a.Position; %gets the positon and size of the color bar
    set(h,'Position',a+[0.05 -0.055 0.008 0.113])% To change size
    if band ==1
        set(h,'XTick',[-0.01 0 0.01])
    elseif band==2
        set(h,...
            'XTick',[-5e-3 0 5e-3])
    elseif band==3
        set(h,...
            'XTick',[-1e-3 0 1e-3])
    end
    export_fig(['figs/rl1_ch', num2str(band),'.png'], '-transparent')
    
    %
    rclean_z = rclean(:,:,band);
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 1000 1000])
    imagesc((rclean_z), clim_log(band,:));
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
    a=get(h);
    a =  a.Position; %gets the positon and size of the color bar
    set(h,'Position',a+[0.05 -0.055 0.008 0.113])% To change size
    if band ==1
        set(h,'XTick',[-0.01 0 0.01])
    elseif band==2
        set(h,...
            'XTick',[-5e-3 0 5e-3])
    elseif band==3
        set(h,...
            'XTick',[-1e-3 0 1e-3])
    end
    export_fig(['figs/rclean_ch', num2str(band),'.png'], '-transparent')
    
    %waitforbuttonpress
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% skip zooms for now %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Plot the west hot spot
% clim_log = [0.008,0.05]; % HS and SARA
% tickLabel = [0.01,0.02,0.03,0.04];
% %clim_log_cl = psf_flux * [0.008 0.05]; % CLEAN
% %tickLabel_cl = psf_flux(band) * [0.01,0.04];
% 
% map_img = cubehelix(2048);
% fontsize=60;
% cen = [1120 400];
% len = 80;
% 
% for band = 1:size(xhs,3)
%     xhs_z = xhs(:,:,band);
%     xhs_z1 = xhs_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
%     f=figure('visible','on');
%     set(gca, 'Color', 'none'); % sets axes background
%     set(f,'PaperUnits','centimeters')
%     set(f,'PaperType','A4');
%     set(f,'PaperOrientation',orient);
%     set(f,'units','pixel','outerposition',[0 0 1250 1250])
%     imagesc((xhs_z1), clim_log);
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
%     set(h,'XTick',tickLabel);
%     a=get(h);
%     a =  a.Position; %gets the positon and size of the color bar
%     set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
%     % Create rectangle
%     annotation(f,'rectangle',...
%     [0.2466 0.106995174129353 0.539200000000001 0.82089552238806],...
%     'Color',[1 1 1],...
%     'LineWidth',4);
%     export_fig(['figs/xhs_z1_ch', num2str(band),'.png'],'-transparent')
%     
%     %
%     xl1_z = xl1(:,:,band);
%     xl1_z1 = xl1_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
%     f=figure('visible','on');
%     set(gca, 'Color', 'none'); % sets axes background
%     set(f,'PaperUnits','centimeters')
%     set(f,'PaperType','A4');
%     set(f,'PaperOrientation',orient);
%     set(f,'units','pixel','outerposition',[0 0 1250 1250])
%     imagesc((xl1_z1), clim_log);
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
%     set(h,'XTick',tickLabel);
%     a=get(h);
%     a =  a.Position; %gets the positon and size of the color bar
%     set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
%     % Create rectangle
%     annotation(f,'rectangle',...
%     [0.2466 0.106995174129353 0.539200000000001 0.82089552238806],...
%     'Color',[1 1 1],...
%     'LineWidth',4);
%     export_fig(['figs/xl1_z1_ch', num2str(band),'.png'], '-transparent')
%     
%     %
%     xclean_z = xclean(:,:,band);
%     xclean_z(xclean_z < 0) = 0;
%     xclean_z1 = xclean_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
%     f=figure('visible','on');
%     set(gca, 'Color', 'none'); % sets axes background
%     set(f,'PaperUnits','centimeters')
%     set(f,'PaperType','A4');
%     set(f,'PaperOrientation',orient);
%     set(f,'units','pixel','outerposition',[0 0 1250 1250])
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
%     export_fig(['figs/xclean_z1_ch', num2str(band),'.png'], '-transparent')
%     
%     %waitforbuttonpress
% end
% 
% %% Plot the inner core
% clim_log = [0.00005  0.0015];
% tickLabel = [1e-4,1e-3];
% %clim_log_cl = psf_flux * [0.00005  0.0015]; % clean
% %tickLabel_cl = round(psf_flux(band) .*[1e-4,1e-3]);
% 
% map_img = cubehelix(2048);
% fontsize=60;
% cen = [N1/2 N2/2+3];
% len = 30;
% 
% for band = 1 : size(xhs,3)
%     xhs_z = xhs(:,:,band);
%     xhs_z1 = xhs_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
%     f=figure('visible','on');
%     set(gca, 'Color', 'none'); % sets axes background
%     set(f,'PaperUnits','centimeters')
%     set(f,'PaperType','A4');
%     set(f,'PaperOrientation',orient);
%     set(f,'units','pixel','outerposition',[0 0 1250 1250])
%     imagesc((xhs_z1), clim_log);
%     colormap(gca, map_img);
%     axis image
%     ax = gca;
%     ax.YAxis.Visible = 'off';
%     ax.XAxis.Visible = 'off';
%     ax.TickLabelInterpreter='latex';
%     h = colorbar;
%     set(h,'Fontsize',fontsize)
%     set(h,'color','white')
%     set(h,'XTick',[1e-4,1e-3]);
%     set(gca,'ColorScale','log');
%     a=get(h);
%     a =  a.Position; %gets the positon and size of the color bar
%     set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
%     % Create ellipse
%     annotation(f,'ellipse',...
%     [0.550090909090909 0.388149939540508 0.0563090909090906 0.0763624982704378],...
%     'Color',[1 1 1],...
%     'LineWidth',5,...
%     'LineStyle','--');
%     % Create rectangle
%     annotation(f,'rectangle',...
%     [0.2466 0.106995174129353 0.539200000000001 0.82089552238806],...
%     'Color',[1 1 1],...
%     'LineWidth',4);
%     export_fig(['figs/xhs_z2_ch', num2str(band),'.png'], '-transparent')
%     
%     %
%     xl1_z = xl1(:,:,band);
%     xl1_z1 = xl1_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
%     f=figure('visible','on');
%     set(gca, 'Color', 'none'); % sets axes background
%     set(f,'PaperUnits','centimeters')
%     set(f,'PaperType','A4');
%     set(f,'PaperOrientation',orient);
%     set(f,'units','pixel','outerposition',[0 0 1250 1250])
%     imagesc((xl1_z1), clim_log);
%     colormap(gca, map_img);
%     axis image
%     ax = gca;
%     ax.YAxis.Visible = 'off';
%     ax.XAxis.Visible = 'off';
%     ax.TickLabelInterpreter='latex';
%     h = colorbar;
%     set(h,'Fontsize',fontsize)
%     set(h,'color','white')
%     set(h,'XTick',[1e-4,1e-3]);
%     set(gca,'ColorScale','log');
%     a=get(h);
%     a =  a.Position; %gets the positon and size of the color bar
%     set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
%     % Create ellipse
%     annotation(f,'ellipse',...
%     [0.550090909090909 0.388149939540508 0.0563090909090906 0.0763624982704378],...
%     'Color',[1 1 1],...
%     'LineWidth',5,...
%     'LineStyle','--');
%     % Create rectangle
%     annotation(f,'rectangle',...
%     [0.2466 0.106995174129353 0.539200000000001 0.82089552238806],...
%     'Color',[1 1 1],...
%     'LineWidth',4);
%     export_fig(['figs/xl1_z2_ch', num2str(band),'.png'], '-transparent')
%     
%     %
%     xclean_z = xclean(:,:,band);
%     xclean_z(xclean_z < 0) = 0;
%     xclean_z1 = xclean_z(cen(1)-len/2:cen(1)+len/2,cen(2)-len/2:cen(2)+len/2);
%     f=figure('visible','on');
%     set(gca, 'Color', 'none'); % sets axes background
%     set(f,'PaperUnits','centimeters')
%     set(f,'PaperType','A4');
%     set(f,'PaperOrientation',orient);
%     set(f,'units','pixel','outerposition',[0 0 1250 1250])
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
%     set(h,'XTick',round(psf_flux(band) .* [1e-4,1e-3],4));
%     %     if band == 1
%     %     set(h,'XTick',[3e-3,5e-2],'XTickLabel',{'3e-3','5e-2'});
%     %     elseif band == 2
%     %     set(h,'XTick',[1e-3,1e-2],'XTickLabel',{'1e-3','1e-2'});
%     %     else
%     %     set(h,'XTick',[1e-4,1e-3],'XTickLabel',{'1e-4','1e-3'});
%     %     end
%     set(gca,'ColorScale','log');
%     a=get(h);
%     a =  a.Position; %gets the positon and size of the color bar
%     set(h,'Position',a+[-0.065 -0.026 0.015 0.055])% To change size
%     % Create ellipse
%     annotation(f,'ellipse',...
%     [0.550090909090909 0.388149939540508 0.0563090909090906 0.0763624982704378],...
%     'Color',[1 1 1],...
%     'LineWidth',5,...
%     'LineStyle','--');
%     % Create rectangle
%     annotation(f,'rectangle',...
%     [0.2466 0.106995174129353 0.539200000000001 0.82089552238806],...
%     'Color',[1 1 1],...
%     'LineWidth',4);
%     export_fig(['figs/xclean_z2_ch', num2str(band),'.png'], '-transparent')
%     
%     %waitforbuttonpress
% end