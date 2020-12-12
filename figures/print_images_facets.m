function print_images_facets(results_path, ncores_data, flag_display, flag_metric)
%%
% Produce the images and metrics reported in the MNRAS paper
% ``A Faceted Prior for Scalable Wideband Imaging: Application to Radio 
% Astronomy''
%
%-------------------------------------------------------------------------%
%%
% Author: P.-A. Thouvenin.
% Last modified: [../../2020]
%-------------------------------------------------------------------------%
%%
% NOTES:
%
% - ...
%-------------------------------------------------------------------------%
%%
% for debugging only
% operatorNorm = 1e4*ones(20,1);
% Qx = [2, 3];
% results_path = '..';
Qx = [1, 2, 3, 4];

%%
% clc; clear all; close all;
format compact;
addpath ../lib/faceted-wavelet-transform/src
addpath ../../CubeHelix
addpath ../../export_fig
savedir = 'figs_facets/';

% Q = Qx.^2;
% doComputeMetric = true; % save metrics in a .txt file
% display = true;

load('ground_truth_spatial_faceting.mat') % .mat file generated with 
                                           % `get_ground_truth.m` -> x0, f,
                                           % operatorNorm
[N(1), N(2), L] = size(x0);
x0 = flipud(x0);             % display image using axis convention from ds9
% overlap_size = 2*floor(N(1)./Qx); % 50% overlap
channels = [1, L];           % channels to be displayed (first and last only)

name_pattern = @(q, sz) fullfile(results_path, ['mnras_faceted_corrected/final_facets/fhs_cw_triangular_N=1024_L=20_Qx=', ...
    num2str(q),'_Qy=',num2str(q),'_Qc=1_ind=1_overlap=',num2str(sz),'_1_1e-05_25.mat']);
mkdir(savedir)

%=========================================================================%
% Reconstruction metrics
%=========================================================================%
% col = @(x) x(:);
norm2D = @(x) squeeze(sqrt(sum(sum(x.^2, 2), 1)));
% norm2 = @(x) squeeze(sqrt(sum(sum(x(:).^2, 2), 1)));
% DR = @(x, res, n) sqrt(prod(N))*squeeze(max(max(x,[],2),[],1))*n./norm2D(res); % per band input
SNR = @(x, x0) 20*log10(norm2D(x0)./norm2D(x - x0));
SNR_log = @(x, x0) SNR(log10(x+eps), log10(x0+eps));

asnr = zeros(numel(Qx), 1);
asnr_log = zeros(numel(Qx), 1);
vsnr = zeros(numel(Qx), 1);
vsnr_log = zeros(numel(Qx), 1);

iteration_number = zeros(numel(Qx), 1);
total_runtime = zeros(numel(Qx), 1);  % total runtime (in s)
total_cpu_time = zeros(numel(Qx), 1); % total CPU time (in s)
aruntime = zeros(numel(Qx), 1);       % average runtime (per iteration, in s)
vruntime = zeros(numel(Qx), 1);       % variance runtime
acpu_time = zeros(numel(Qx), 1);      % average cpu time (per iter., in s)
vcpu_time = zeros(numel(Qx), 1);      % variance cpu time
atime_facet = zeros(numel(Qx), 1);    % average time (facet, per iter., in s)
atime_data = zeros(numel(Qx), 1);     % variance
vtime_facet = zeros(numel(Qx), 1);    % average time (data, per iter., in s)
vtime_data = zeros(numel(Qx), 1);     % variance

%=========================================================================%
% Plot parameters
%=========================================================================%
clim = [1e-5, 1;  % image
    1e-4, 1;            % error image
    -3.5e-6, 3.5e-6]; % residual images % e-4 before

clim(:,:,2) = [5e-5, 1; % image
    1e-3, 1;                % error image
    -3.5e-6, 3.5e-6];     % residual images

fontsize=25;
map_img = cubehelix(256);

%% Display ground truth image (spatial faceting)
if flag_display
    for l = 1:2
        % ground truth
        display_image(x0(:,:,channels(l)), clim(1,:,l), map_img, fontsize, true);
        export_fig(strcat(savedir,'x', num2str(l),'.pdf'), '-transparent','-q101')
        close
    end
end

%% Load images %! first and last channel only for SARA, overall norm otherwise?
for k = 1:numel(Qx)
    if Qx(k) > 1
        % faceted HyperSARA
        if flag_metric || (flag_display && Qx(k) == 4)
            overlap = floor(N(1)/Qx(k));
            fileName = name_pattern(Qx(k), overlap);
            load(fileName, 'xsol', 'res', 't_facet', 't_data', 'end_iter')
            x = flipud(xsol);
        end

        % snr
        if flag_metric
            a = SNR(x, x0);
            asnr(k) = mean(a);
            vsnr(k) = var(a);
            a = SNR_log(x, x0);
            asnr_log(k) = mean(a);
            vsnr_log(k) = var(a);

            % timing
            iteration_number(k) = sum((end_iter > 0));
            %
            total_runtime(k) = sum(end_iter(end_iter > 0));
            aruntime(k) = total_runtime(k)/iteration_number(k); % average runtime per iteration
            vruntime(k) = var(end_iter(end_iter > 0));
            %
            atime_facet(k) = mean(t_facet(t_facet > 0));
            vtime_facet(k) = var(t_facet(t_facet > 0));
            atime_data(k) = mean(t_data(t_data > 0));
            vtime_data(k) = var(t_data(t_data > 0));
            %
            a_ = Qx(k)^2*t_facet(t_facet > 0) + ncores_data*t_data(t_data > 0);
            total_cpu_time(k) = sum(a_);
            acpu_time(k) = Qx(k)^2*atime_facet(k) + ncores_data*atime_data(k);
            vcpu_time(k) = var(a_);
        end
        
        if flag_display && Qx(k) == 4
            res = flipud(res(:,:,[1,end]))./reshape(operatorNorm([1,end]), [1, 1, 2]);
            x = x(:,:,[1,end]); % keep last two channels
            for l = 1:2
                % images
                display_image(x(:,:,l), clim(1,:,l), map_img, fontsize, true);
                export_fig(strcat(savedir,'x_fhs', num2str(l),'.pdf'), '-transparent','-q101')
                close

                % residual images
                display_image(res(:,:,l), clim(3,:,l), map_img, fontsize, false);
                export_fig(strcat(savedir,'res_fhs', num2str(l),'.pdf'), '-transparent','-q101')
                close
            end
        end
    else
        % HyperSARA results
        if flag_metric || flag_display
            fileName = fullfile(results_path, 'mnras_faceted_corrected/final_hypersara/fhs_hypersara_triangular_N=1024_L=20_Qx=1_Qy=1_Qc=1_ind=1_overlap=256_1_1e-05_25.mat');
            load(fileName, 'xsol', 'res', 't_l21', 't_nuclear', 't_master', 't_data', 'end_iter')
            x = flipud(xsol);
        end

        % snr
        if flag_metric
            a = SNR(x, x0);
            asnr(k) = mean(a);
            vsnr(k) = var(a);
            a = SNR_log(x, x0);
            asnr_log(k) = mean(a);
            vsnr_log(k) = var(a);

            % timing
            a_ = t_l21(t_l21 > 0) + t_nuclear(t_nuclear > 0) + t_master(t_master > 0);
            b_ = t_data(t_data > 0);
            iteration_number(k) = sum((end_iter > 0));
            %
            atime_facet(k) = mean(a_);
            vtime_facet(k) = var(a_);
            atime_data(k) = mean(b_);
            vtime_data(k) = var(b_);
            %
            total_runtime(k) = sum(end_iter(end_iter > 0));
            aruntime(k) = total_runtime(k)/iteration_number(k); % average runtime per iteration
            vruntime(k) = var(end_iter(end_iter > 0));
            %
            a_ = a_ + ncores_data*b_;
            total_cpu_time(k) = sum(a_);
            acpu_time(k) = atime_facet(k) + ncores_data*atime_data(k); % average cpu time per iteration
            vcpu_time(k) = var(a_);
        end
        
        % display image
        if flag_display
            res = flipud(res(:,:,[1,end]))./reshape(operatorNorm([1,end]), [1, 1, 2]);
            x = x(:,:,[1,end]); % keep last two channels
            for l = 1:2
                % images
                display_image(x(:,:,l), clim(1,:,l), map_img, fontsize, true);
                export_fig(strcat(savedir,'x_hs', num2str(l),'.pdf'), '-transparent','-q101')
                close

                % residual images
                display_image(res(:,:,l), clim(3,:,l), map_img, fontsize, false);
                export_fig(strcat(savedir,'res_hs', num2str(l),'.pdf'), '-transparent','-q101')
                close
            end
        end
    end
end

%% Display results (table)
if flag_metric
    for k = 1:numel(Qx)
       fprintf("Qx = %i, asnr = %2.2f, std_snr = %1.2e, asnr_log = %2.2f, std_snr_log = %1.2e, iteration_number = %i \n", ...
           Qx(k), asnr(k), std(vsnr(k)), asnr_log(k), std(vsnr_log(k)), iteration_number(k))
       fprintf(" aruntime (s) = %.2f, vruntime (s) = %1.2e, acpu_time (s) = %.2f, vcpu_time (s) = %1.2e \n", ...
           aruntime(k), std(vruntime(k)), acpu_time(k), std(vcpu_time(k)));
       fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
           total_runtime(k)/3600, total_cpu_time(k)/3600)
    end

%% Saving results
    save('results_facets.mat', '-v7.3', 'asnr', 'vsnr', 'asnr_log', ...
        'vsnr_log', 'aruntime', 'vruntime', 'acpu_time', 'vcpu_time', ...
        'atime_facet', 'vtime_facet', 'atime_data', 'vtime_data', ...
        'iteration_number', 'total_runtime', 'total_cpu_time');
end

%%
% %% Faceted HyperSARA (results with mutliple number of facets and overlap)
% fileNameFunction = @(Q, d) strcat('new_results/results_hyperSARA_spmd4_cst_weighted_img=1024_tot=20_Qx=', ...
%     num2str(Q), '_Qy=', num2str(Q), '_d=', num2str(d), '_Qc=20_gamma=1e-05', '.mat');
% 
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
%         fileName_tmp = fileNameFunction(nfacets(n), d{n}(k));
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
%         for l = 1:2
%             % images
%             f=figure('visible','on');
%             set(gca, 'Color', 'none'); % sets axes background
%             set(f,'PaperUnits','centimeters')
%             set(f,'PaperType','A4');
%             set(f,'PaperOrientation',orient);
%             set(f,'units','pixel','outerposition',[0 0 600 600])
%             imagesc(log10(x_facet(:,:,l)), clim_log(1,:,l));
%             colormap(gca, map_img);
%             axis image
%             ax = gca;
%             ax.YAxis.Visible = 'off';
%             ax.XAxis.Visible = 'off';
%             h = colorbar;
%             set(h,'Fontsize',fontsize);
%             export_fig(strcat(savedir,'x_fhs', num2str(l), '_Q', num2str(nfacets(n)),...
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
%             imagesc(res_facet(:,:,l), clim_log(3,:,l));
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
%             export_fig(strcat(savedir,'res_fhs', num2str(l), '_Q', num2str(nfacets(n)),...
%                 '_d', num2str(d{n}(k)),'.pdf'), '-transparent','-q101')
%             
%             close 
%         end
%     end
% end