clc; clear all; close all;
format compact;

addpath lib/sdwt2
addpath ../CubeHelix
addpath ../export_fig
mkdir figs

%%
% Produce the images reported in the reesponse to the 2020 MNRAS paper.
%
% Results with the triangular apodization window for the faceted approach.
%-------------------------------------------------------------------------%
%%
% Author: P.-A. Thouvenin. [21/09/2020]
% Last modified: [21/09/2020]
%-------------------------------------------------------------------------%
%%
% NOTES:
%
% - ...
%-------------------------------------------------------------------------%
%%
savedir = 'figs/M31/';

%% Load ground truth image

% parameters chaging from one run to another
% N = [1024,1024];
% Qx = 4;
% Qy = 4;
% d = {256}; % check overlap (to be included in the bale of the file)
% Q = Qx*Qy;
% algo_version = 'cst_weighted';
% 
% % parameters kept constant over all simulations
% window_type = 'triangular';
% Qc = 1;
% nchannels = 20;
% p = 1;
% snr = 60;
% 
% results_path = 'results/W28_1024_m39_1';

% M31
N = [256,256];
Qx = 2;
Qy = 1;
d = {128}; % check overlap (to be included in the bale of the file)
Q = Qx*Qy;
algo_version = 'cst_weighted';

% parameters kept constant over all simulations
window_type = 'triangular';
Nx = 256;
Qc = 1;
nchannels = 15;
p = 1;
snr = 60;

results_path = 'results/M31';

data_path = 'data';

%%
% TO BE MODIFIED LATER ON
% fileNameFunction = @(Qx, Qy, Qc) strcat('facetHyperSARA_cst_weighted_triangular_N=256_L=15_p=1_Qx=', ...
%     num2str(Qx),'_Qy=',num2str(Qy),'_Qc=', num2str(Qc), '_snr=60.mat');
% load(fullfile(results_path,fileNameFunction(Qx,Qy,Qc)), 'xsol', 'param');

% M31
load(fullfile(results_path,'facetHyperSARA_cst_weighted_triangular_N=256_L=15_p=1_Qx=2_Qy=1_Qc=1_snr=60.mat'), 'xsol', 'param');
load(fullfile(results_path,'facetHyperSARA_c_w_0_0.01_30.mat'), 'res');
x = fitsread(fullfile(data_path,'M31_L20.fits')).';

% W28_m39_1
% load(fullfile(results_path,'fhs_cst_weighted_triangular_N=1024_L=20_p=0.5_Qx=4_Qy=4_Qc=1_overlap=128_snr=60.mat'), 'xsol', 'param');
% load(fullfile(results_path,'facetHyperSARA_c_w_W28_1024_m39_0_0.0056566_0.mat'), 'res');
% x = fitsread(fullfile(data_path,'W28_512_m39_L20.fits')).';


% results_name_function = @(Qx, Qy, overlap) strcat('fhs_', algo_version,'_',window_type,'_N=',num2str(Nx), ...
%     '_L=',num2str(nchannels),'_p=',num2str(p), ...
%     '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc), ...
%     '_overlap=', num2str(overlap_size), ...
%     '_snr=', num2str(input_snr),'.mat');
% 
% load(fullfile(results_path,'facetHyperSARA_c_w_0_0.01_1.mat'), 'res');

Nx = sqrt(size(x, 1));
Ny = Nx;
nChannels = size(x, 2);
x = reshape(x, [Ny, Nx, nChannels]);
x = flipud(x(:,:,[1, nchannels]));
xsol = flipud(xsol(:,:,[1, nchannels]));
res = flipud(res(:,:,[1, nchannels])/param.nu2);

%% Diplay results (images, just first and last band)
% ground truth, hyperSARA, SARA

% plot residual image faceted HyperSARA for multiple values of Q, and
% multiple values of overlap (..., to be specified beforehand)

%=========================================================================%
% Plot parameters
%=========================================================================%

clim_log = [-5, ceil(log10(max(max(xsol(:,:,1)))));    % image
    -0.7*max(max(res(:,:,1))), 0.7*max(max(res(:,:,1)))]; % residual images

clim_log(:,:,2) = [-4, ceil(log10(max(max(xsol(:,:,2)))));   % image
    -0.7*max(max(res(:,:,2))), 0.7*max(max(res(:,:,2)))]; % residual images

% M31
% clim_log = [-5, 1;    % image
%     -3.5e-2, 3.5e-2]; % residual images
% 
% clim_log(:,:,2) = [-4, 1;   % image
%     -3.5e-2, 3.5e-2];       % residual images

% clim_log = [-5, 0];    % image
%     -3.5e-4, 3.5e-4]; % residual images
% 
% clim_log(:,:,2) = [-4, 0;   % image
%     -3.5e-4, 3.5e-4];       % residual images

fontsize=25;
map_img = cubehelix(256); % colormap for the images
map = cubehelix(256);     % colormap for the residual image
linewidth = 2.3;

mkdir(savedir)

%%
 % non-overlapping tessellation
 rg_y = domain_decomposition(Qy, N(1));
 rg_x = domain_decomposition(Qx, N(2));
 I = zeros(Q, 2);
 dims = zeros(Q, 2);
 for qx = 1:Qx
     for qy = 1:Qy
         q = (qx-1)*Qy+qy;
         I(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1];
         dims(q, :) = [rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
     end
 end
 
 % redundant facets
 rg_yo = domain_decomposition_overlap2(Qy, N(1), d{1}(1));
 rg_xo = domain_decomposition_overlap2(Qx, N(2), d{1}(1));
 Io = zeros(Q, 2);
 dims_o = zeros(Q, 2);
 overlap = zeros(Q, 2);
 for qx = 1:Qx
     for qy = 1:Qy
         q = (qx-1)*Qy+qy;
         Io(q, :) = [rg_yo(qy, 1)-1, rg_xo(qx, 1)-1];
         dims_o(q, :) = [rg_yo(qy,2)-rg_yo(qy,1)+1, rg_xo(qx,2)-rg_xo(qx,1)+1];
         overlap(q, :) =  dims_o(q,:) - dims(q,:);
     end
 end
 ds = dims_o - dims;
 c = floor(Qx/2)*Qy + floor(Qy/2) + 1; % index of the central facet

 %%
for l = 1:2
    % ground truth
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 600 600])
    imagesc(log10(x(:,:,l)), clim_log(1,:,l));
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    export_fig(strcat(savedir,'x', num2str(l),'.pdf'), '-transparent','-q101')
    close
    
    %--
    % reconstructed
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 600 600])
    imagesc(log10(xsol(:,:,l)), clim_log(1,:,l));
    colormap(gca, map_img);
    axis image
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    h = colorbar;
    set(h,'Fontsize',fontsize)
    export_fig(strcat(savedir,'x_fhs', num2str(l),'.pdf'), '-transparent','-q101')
    close
    
    %--
    % residual images
    f=figure('visible','on');
    set(gca, 'Color', 'none'); % sets axes background
    set(f,'PaperUnits','centimeters')
    set(f,'PaperType','A4');
    set(f,'PaperOrientation',orient);
    set(f,'units','pixel','outerposition',[0 0 600 600])
    imagesc(res(:,:,l), clim_log(2,:,l));
    colormap(gca, map);
    axis image
    h = colorbar;
    set(h,'Fontsize',fontsize)
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'off';
    % horizontal lines
    for qx = 1:Qx-1
        q = (qx-1)*Qy + qx;
        xx = [1, N(1)];
        yy = (I(q, 1) + dims(q, 1))*ones(1, 2);
        line(xx,yy,'Color','r','LineStyle','--','LineWidth',linewidth)
    end
    % vertical lines
    for qx = 1:Qx-1
        q = (qx-1)*Qy + qx;
        yy = [1, N(2)];
        xx = (I(q, 2) + dims(q, 2))*ones(1, 2);
        line(xx,yy,'Color','r','LineStyle','--','LineWidth',linewidth)
    end
    % central facet
    rectangle('Position', [I(c, 2)-ds(c, 2), I(c, 1)-ds(c, 1), ...
        dims_o(c, 2), dims_o(c, 1)], ...
        'LineWidth',linewidth,'LineStyle','-','EdgeColor','r')
    export_fig(strcat(savedir,'res_fhs', num2str(l), '_Q', num2str(Q),...
        '_d', num2str(d{1}(1)),'.pdf'), '-transparent','-q101')
    
    close
end

%% Faceted HyperSARA (results with mutliple number of facets and overlap)
% fileNameFunction = @(Q, d) strcat('new_results/results_hyperSARA_spmd4_cst_weighted_img=1024_tot=20_Qx=', ...
%     num2str(Q), '_Qy=', num2str(Q), '_d=', num2str(d), '_Qc=20_gamma=1e-05', '.mat');

% nfacets = [2, 3, 4];
% d = cell(size(nfacets));
% d{1} = [128; 256; 512];
% d{2} = [114; 171; 342];
% d{3} = [16, 64, 128, 256];
% 
% % facet hyperSARA
% for n = 1:numel(nfacets)
%     
%     Qx = nfacets(n);
%     Qy = Qx;
%     Q = Qx*Qy;
%     for k = 1:numel(d{n})
%         
%         % images to be kept in memory (for the next section)
%         % fileName_tmp = fileNameFunction(nfacets(n), d{n}(k));
%         
%         load(fileName_tmp, 'xsol', 'res')
%         x_facet = flipud(xsol);
%         res_facet = flipud(res(:,:,[1,end]))./reshape(norm_Phi, [1, 1, 2]);
%         err_facet = abs(x0(:,:,[1,end]) - x_facet(:,:,[1,end]));
%         x_facet = x_facet(:,:,[1, end]);
%         
%         % non-overlapping tessellation
%         rg_y = domain_decomposition(Qy, N(1));
%         rg_x = domain_decomposition(Qx, N(2));
%         I = zeros(Q, 2);
%         dims = zeros(Q, 2);
%         for qx = 1:Qx
%             for qy = 1:Qy
%                 q = (qx-1)*Qy+qy;
%                 I(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1];
%                 dims(q, :) = [rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
%             end
%         end
% 
%         % redundant facets
%         rg_yo = domain_decomposition_overlap2(Qy, N(1), d{n}(k));
%         rg_xo = domain_decomposition_overlap2(Qx, N(2), d{n}(k));
%         Io = zeros(Q, 2);
%         dims_o = zeros(Q, 2);
%         overlap = zeros(Q, 2);
%         for qx = 1:Qx
%             for qy = 1:Qy
%                 q = (qx-1)*Qy+qy;
%                 Io(q, :) = [rg_yo(qy, 1)-1, rg_xo(qx, 1)-1];
%                 dims_o(q, :) = [rg_yo(qy,2)-rg_yo(qy,1)+1, rg_xo(qx,2)-rg_xo(qx,1)+1];
%                 overlap(q, :) =  dims_o(q,:) - dims(q,:);%max(dims_overlap_ref(q, :), dims_o(q, :)) - dims(q,:); % issue here! max(dims_overlap{q}, [], 1)
%             end
%         end
%         ds = dims_o - dims;
%         c = floor(Qx/2)*Qy + floor(Qy/2) + 1; % index of the central facet
% 
%         for i = 1:2
%             % images
%             f=figure('visible','on');
%             set(gca, 'Color', 'none'); % sets axes background
%             set(f,'PaperUnits','centimeters')
%             set(f,'PaperType','A4');
%             set(f,'PaperOrientation',orient);
%             set(f,'units','pixel','outerposition',[0 0 600 600])
%             imagesc(log10(x_facet(:,:,i)), clim_log(1,:,i));
%             colormap(gca, map_img);
%             axis image
%             ax = gca;
%             ax.YAxis.Visible = 'off';
%             ax.XAxis.Visible = 'off';
%             h = colorbar;
%             set(h,'Fontsize',fontsize);
%             export_fig(strcat(savedir,'x_fhs', num2str(i), '_Q', num2str(nfacets(n)),...
%                 '_d', num2str(d{n}(k)),'.pdf'), '-transparent','-q101')
%             close
%             %--
%             % residual images
%             f=figure('visible','on');
%             set(gca, 'Color', 'none'); % sets axes background
%             set(f,'PaperUnits','centimeters')
%             set(f,'PaperType','A4');
%             set(f,'PaperOrientation',orient);
%             set(f,'units','pixel','outerposition',[0 0 600 600])
%             imagesc(res_facet(:,:,i), clim_log(3,:,i));
%             colormap(gca, map);
%             axis image
%             h = colorbar;
%             set(h,'Fontsize',fontsize)
%             ax = gca;
%             ax.YAxis.Visible = 'off';
%             ax.XAxis.Visible = 'off';
%             % horizontal lines
%             for qx = 1:Qx-1
%                 q = (qx-1)*Qy + qx;
%                 xx = [1, N(1)];
%                 yy = (I(q, 1) + dims(q, 1))*ones(1, 2);
%                 line(xx,yy,'Color','r','LineStyle','--','LineWidth',linewidth)
%             end
%             % vertical lines
%             for qx = 1:Qx-1
%                 q = (qx-1)*Qy + qx;
%                 yy = [1, N(2)];
%                 xx = (I(q, 2) + dims(q, 2))*ones(1, 2);
%                 line(xx,yy,'Color','r','LineStyle','--','LineWidth',linewidth)
%             end
%             % central facet
%             rectangle('Position', [I(c, 2)-ds(c, 2), I(c, 1)-ds(c, 1), ...
%                 dims_o(c, 2), dims_o(c, 1)], ...
%                 'LineWidth',linewidth,'LineStyle','-','EdgeColor','r')
%             export_fig(strcat(savedir,'res_fhs', num2str(i), '_Q', num2str(nfacets(n)),...
%                 '_d', num2str(d{n}(k)),'.pdf'), '-transparent','-q101')
%             
%             close 
%         end
%     end
% end
