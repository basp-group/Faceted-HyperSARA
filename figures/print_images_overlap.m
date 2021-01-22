function print_images_overlap(results_path, ncores_data)
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

asnr = zeros(numel(overlap_size), 1);
asnr_log = zeros(numel(overlap_size), 1);
vsnr = zeros(numel(overlap_size), 1);
vsnr_log = zeros(numel(overlap_size), 1);

iteration_number = zeros(numel(overlap_size), 1);
total_runtime = zeros(numel(overlap_size), 1);  % total runtime (in s)
total_cpu_time = zeros(numel(overlap_size), 1); % total CPU time (in s)
aruntime = zeros(numel(overlap_size), 1);       % average runtime (per iteration, in s)
vruntime = zeros(numel(overlap_size), 1);       % variance runtime
acpu_time = zeros(numel(overlap_size), 1);      % average cpu time (per iter., in s)
vcpu_time = zeros(numel(overlap_size), 1);      % variance cpu time
atime_facet = zeros(numel(overlap_size), 1);    % average time (facet, per iter., in s)
atime_data = zeros(numel(overlap_size), 1);     % variance
vtime_facet = zeros(numel(overlap_size), 1);    % average time (data, per iter., in s)
vtime_data = zeros(numel(overlap_size), 1);     % variance

%=========================================================================%
% Plot parameters
%=========================================================================%
% clim_log = [1e-5, 1;  % image
%     1e-4, 1;            % error image
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
    
    load(fileName, 'xsol', 't_facet', 't_data', 'end_iter') % 'res'
    
    % snr
    x = flipud(xsol);
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
    vtime_facet(k) = (Q^2)*var(t_facet(t_facet > 0)); % t_facet already averaged over facet cores -> multiplicative factor to retrieve the true variance
    atime_data(k) = mean(t_data(t_data > 0));
    vtime_data(k) = (ncores_data^2)*var(t_data(t_data > 0));
    %
    a_ = Q*t_facet(t_facet > 0) + ncores_data*t_data(t_data > 0); % average cpu time per iteration
    total_cpu_time(k) = sum(a_);
    acpu_time(k) = total_cpu_time(k)/iteration_number(k);
    vcpu_time(k) = var(a_);
    
%     res = flipud(res(:,:,[1,end]))./reshape(operatorNorm, [1, 1, L]);
%     x = x(:,:,[1,end]); % keep last two channels
%     
%     for l = 1:2
%         % images
%         display_image(x(:,:,l), clim_log(1,:,l), map_img, fontsize, true);
%         export_fig(strcat(savedir,'x_fhs', num2str(l),'.pdf'), '-transparent','-q101')
%         close
%         
%         % residual images
%         display_image(res(:,:,l), clim_log(3,:,l), map_img, fontsize, false);
%         export_fig(strcat(savedir,'res_hs', num2str(l),'.pdf'), '-transparent','-q101')
%         close
%     end
end

%% Display results (table)
for k = 1:numel(overlap_size)
    fprintf("overlap = %i, asnr = %2.2f, std_snr = %1.2e, asnr_log = %2.2f, std_snr_log = %1.2e, iteration_number = %i \n", ...
       overlap_size(k), asnr(k), sqrt(vsnr(k)), asnr_log(k), sqrt(vsnr_log(k)), iteration_number(k))
   fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
       aruntime(k), sqrt(vruntime(k)), acpu_time(k), sqrt(vcpu_time(k)));
   fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
       total_runtime(k)/3600, total_cpu_time(k)/3600)
end

%% Saving results
save('results_overlap.mat', '-v7.3', 'asnr', 'vsnr', 'asnr_log', ...
    'vsnr_log', 'aruntime', 'vruntime', 'acpu_time', 'vcpu_time', ...
    'atime_facet', 'vtime_facet', 'atime_data', 'vtime_data', ...
    'iteration_number', 'total_runtime', 'total_cpu_time');
