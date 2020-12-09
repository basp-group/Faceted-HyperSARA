function print_images_overlap(results_path)
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
%% Debugging
% results_path = '..';
% % for debugging only
% operatorNorm = 1e4*ones(L,1);
% overlap_size = [20]; % 

%%
format compact;
addpath ../lib/faceted-wavelet-transform/src
addpath ../../CubeHelix
addpath ../../export_fig

overlap_size = [0; 20; 64; 128];
Qx = 4;
Qy = 4;
Q = Qx*Qy;
% doComputeMetric = true; % save metrics in a .txt file
% display = true;

load('ground_truth_spatial_faceting.mat') % .mat file generated with 
                                           % `get_ground_truth.m` -> x0, f,
                                           % operatorNorm                                           
[N(1), N(2), L] = size(x0);
x0 = flipud(x0);             % display image using axis convention from ds9

file_pattern = @(sz) fullfile(results_path,['mnras_faceted_corrected/final_overlap/fhs_cw_triangular_N=1024_L=20_Qx=4_Qy=4_Qc=1_ind=1_overlap=', ...
    num2str(sz), '_1_1e-05_25.mat']);
file_pattern2 = @(sz) fullfile(results_path,['mnras_faceted_corrected/final_overlap/fhs_no_triangular_N=1024_L=20_Qx=4_Qy=4_Qc=1_ind=1_overlap=', ...
    num2str(sz), '_1_1e-05_25.mat']);
% mkdir(savedir)

%=========================================================================%
% Reconstruction metrics
%=========================================================================%
col = @(x) x(:);
norm2D = @(x) squeeze(sqrt(sum(sum(x.^2, 2), 1)));
norm2 = @(x) squeeze(sqrt(sum(sum(x(:).^2, 2), 1)));
DR = @(x, res, n) sqrt(prod(N))*squeeze(max(max(x,[],2),[],1))*n./norm2D(res); % per band input
SNR = @(x, x0) 20*log10(norm2D(x0)./norm2D(x - x0));
SNR_log = @(x, x0) SNR(log10(x+eps), log10(x0+eps));

asnr = zeros(numel(overlap_size));
asnr_log = zeros(numel(overlap_size));
vsnr = zeros(numel(overlap_size));
vsnr_log = zeros(numel(overlap_size));
runtime = zeros(numel(overlap_size));
cpu_time = zeros(numel(overlap_size));
atime_facet = zeros(numel(overlap_size), 1);
atime_data = zeros(numel(overlap_size), 1);
vtime_facet = zeros(numel(overlap_size), 1);
vtime_data = zeros(numel(overlap_size), 1);

%=========================================================================%
% Plot parameters
%=========================================================================%
% clim_log = [-5, 0;  % image
%     -4, 0;            % error image
%     -3.5e-4, 3.5e-4]; % residual images
% 
% clim_log(:,:,2) = [-4, 0; % image
%     -3, 0;                % error image
%     -3.5e-4, 3.5e-4];     % residual images
% 
% fontsize=25;
% map_img = cubehelix(256);
% linewidth = 2.3;

%% Load images
for k = 1:numel(overlap_size)
    if overlap_size(k) == 0
        fileName = file_pattern2(overlap_size(k));
    else
        fileName = file_pattern(overlap_size(k));
    end
    
    load(fileName, 'xsol', 't_facet', 't_data', 'end_iter') % 'res', 'param'
    % operator_norm = param.Anorm;
    x = flipud(xsol);
    
    a = SNR(x, x0);
    asnr(k) = mean(a);
    vsnr(k) = var(a);
    a = SNR_log(x, x0);
    asnr_log(k) = mean(a);
    vsnr_log(k) = var(a);
    
    runtime(k) = mean(end_iter(end_iter > 0)); % average runtime per iteration
    atime_facet(k) = mean(t_facet(t_facet > 0));
    atime_data(k) = mean(t_data(t_data > 0));
    vtime_facet(k) = var(t_facet(t_facet > 0));
    vtime_data(k) = var(t_data(t_data > 0));
    cpu_time(k) = Q*atime_facet(k) + L*atime_data(k); % average cpu time per iteration
    
%     res = flipud(res(:,:,[1,end]))./reshape(operatorNorm, [1, 1, L]);
%     x = x(:,:,[1,end]); % keep last two channels
%     
%     for l = 1:2
%         % images
%         display_image(log10(x(:,:,l)), clim_log(1,:,l), map_img, fontsize);
%         export_fig(strcat(savedir,'x_fhs', num2str(l),'.pdf'), '-transparent','-q101')
%         close
%         
%         % residual images
%         display_image(res(:,:,l), clim_log(3,:,l), map_img, fontsize);
%         export_fig(strcat(savedir,'res_hs', num2str(l),'.pdf'), '-transparent','-q101')
%         close
%     end
end

%% Display results (table)
for k = 1:numel(overlap_size)
   fprintf("overlap = %i, asnr = %2.2f, vsnr = %1.2e, asnr_log = %2.2f, vsnr_log = %1.2e, runtime = %.2f \n", ...
       overlap_size(k), asnr(k), vsnr(k), asnr_log(k), vsnr_log(k), runtime(k))
end

%% Saving results
save('results_overlap.mat', '-v7.3', 'asnr', 'vsnr', 'asnr_log', ...
    'vsnr_log', 'runtime', 'cpu_time', 'atime_facet', 'atime_data', ...
    'vtime_facet', 'vtime_data')
