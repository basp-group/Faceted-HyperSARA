function print_images_spectral(results_path, ncores_data, flag_display, flag_metric)
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
% - display first channel from the first block, last channel of the last
% block for Qc = 10. Display only HyperSARA results (not SARA).
% - need to retrive ground truth images first (run get_ground_truth.m)
% - need to renormalize residual images with the operator norm computed for
% each channel independently (?)
%
% TODO: load operator norm, compute snr (print in a .txt file), compute
% timing, relate to number of CPUs
%
%-------------------------------------------------------------------------%
%%
% debugging
% results_path = '..';
% operatorNorm = 1e4*ones(100,1);
% flag_metric = true;
% flag_display = true;

%%
% clc; clear all; close all;
format compact;
addpath ../lib/faceted-wavelet-transform/src
addpath ../../CubeHelix
addpath ../../export_fig
savedir = 'figs_spectral/';
% Qx = 1;
% Qy = 1;
% Q = Qx*Qy;
Qc = [1, 2, 3, 5, 7, 10, 16];

load('ground_truth_spectral_faceting.mat') % .mat file generated with 
                                           % `get_ground_truth.m` -> x0, f,
                                           % operatorNorm
[N(1), N(2), L] = size(x0);
x0 = flipud(x0);             % display image using axis convention from ds9
chans = [1, 100];

name_pattern = @(Qc, ind) fullfile(results_path, ['mnras_faceted_corrected/final_spectral/fhs_hypersara_triangular_N=', ...
    num2str(N(1)),'_L=',num2str(L),'_Qx=1_Qy=1_Qc=', ...
    num2str(Qc),'_ind=',num2str(ind),'_overlap=256_', ...
    num2str(ind),'_0.001_10.mat']);
mkdir(savedir)

%=========================================================================%
% Reconstruction metrics
%=========================================================================%
% col = @(x) x(:);
norm2D = @(x) squeeze(sqrt(sum(sum(x.^2, 2), 1)));
% norm2 = @(x) squeeze(sqrt(sum(sum(x(:).^2, 2), 1)));
% DR = @(x, res, n) sqrt(prod(N))*squeeze(max(max(x,[],2),[],1))*n./norm2D(res); % per band input
SNR = @(x, x0) 20*log10(norm2D(x0)./norm2D(x - x0));
% SNR_log = @(x, x0) SNR(log10(x+eps), log10(x0+eps));
SNR_log = @(x, x0) 20*log10(norm2D(log10(x0+eps))./norm2D(log10(x+eps)-log10(x0+eps)));

asnr = zeros(numel(Qc), 1);
asnr_log = zeros(numel(Qc), 1);
vsnr = zeros(numel(Qc), 1);
vsnr_log = zeros(numel(Qc), 1);

total_runtime = zeros(numel(Qc), 1);  % assuming all sub-cubes are running in parallel
total_cpu_time = zeros(numel(Qc), 1); % 
aruntime = zeros(numel(Qc), 1);
vruntime = zeros(numel(Qc), 1);
acpu_time = zeros(numel(Qc), 1);
vcpu_time = zeros(numel(Qc), 1);
atime_l21 = zeros(numel(Qc), 1); % average time (per iteration)
atime_nuclear = zeros(numel(Qc), 1);
atime_data = zeros(numel(Qc), 1);
atime_master = zeros(numel(Qc), 1);
vtime_l21 = zeros(numel(Qc), 1); % variance
vtime_nuclear = zeros(numel(Qc), 1);
vtime_data = zeros(numel(Qc), 1);
vtime_master = zeros(numel(Qc), 1);
iteration_number = zeros(numel(Qc), 1);

%=========================================================================%
% Plot parameters
%=========================================================================%
clim_log = [1e-5, 1;  % image
    1e-4, 1;            % error image
    -3.5e-6, 3.5e-6]; % residual images % e-4 before

clim_log(:,:,2) = [5e-5, 1; % image
    1e-3, 1;                % error image
    -3.5e-6, 3.5e-6];     % residual images

fontsize = 25;
map_img = cubehelix(256);

%% Display ground truth image (spectral faceting)
if flag_display
    for l = 1:2
        % ground truth
        display_image(x0(:,:,chans(l)), clim_log(1,:,l), map_img, fontsize, true);
        export_fig(strcat(savedir,'x', num2str(l),'_sub.pdf'), '-transparent','-q101')
        close
    end
end

%% Results (HyperSARA, spectrally faceted HyperSARA)
for k = 1:numel(Qc)
    x = zeros(N(1), N(2), L);
    res_ = zeros(N(1), N(2), L);
    channels = split_range_interleaved(Qc(k), L);
    
    sum_runtime_sqr = 0;
    sum_l21_sqr = 0;
    sum_nuclear_sqr = 0;
    sum_master_sqr = 0;
    sum_data_sqr = 0;
    sum_cpu_sqr = 0;
    for ind = 1:Qc(k)
        fileName = name_pattern(Qc(k), ind);
        % ncpus =
        load(fileName, 'xsol', 'res', 't_l21', 't_nuclear', 't_master', 't_data', 'end_iter')
        % operator_norm = param.nu2;
        x(:,:,channels{ind}) = flipud(xsol);
        res_(:,:,channels{ind}) = flipud(res); 
        
        if flag_metric
            % compute mean and variance over all the files? just 1 for now
            aruntime(k) = aruntime(k) + sum(end_iter(end_iter > 0)); % average runtime per iteration, over all sub-problems
            sum_runtime_sqr = sum_runtime_sqr + sum(end_iter(end_iter > 0).^2);

            atime_l21(k) = atime_l21(k) + sum(t_l21(t_l21 > 0));
            atime_nuclear(k) = atime_nuclear(k) + sum(t_nuclear(t_nuclear > 0));
            atime_master(k) = atime_master(k) + sum(t_master(t_master > 0));
            atime_data(k) = atime_data(k) + sum(t_data(t_data > 0));
            sum_l21_sqr = sum_l21_sqr + sum(t_l21(t_l21 > 0).^2);
            sum_nuclear_sqr = sum_nuclear_sqr + sum(t_nuclear(t_nuclear > 0).^2);
            sum_master_sqr = sum_master_sqr + sum(t_master(t_master > 0).^2);
            sum_data_sqr = sum_data_sqr + sum(t_data(t_data > 0).^2);

            % average number of iterations over all sub-problems
            iteration_number(k) = iteration_number(k) + sum(end_iter > 0);
            total_runtime(k) = max(total_runtime(k), sum(end_iter(end_iter > 0)));
            total_cpu_time(k) = total_cpu_time(k) + sum(t_l21(t_l21 > 0)) + ...
                sum(t_nuclear(t_nuclear > 0)) + sum(t_master(t_master > 0)) ...
                + ncores_data*sum(t_data(t_data > 0));
            sum_cpu_sqr = sum_cpu_sqr + sum((t_l21(t_l21 > 0) + t_nuclear(t_nuclear > 0) + t_master(t_master > 0) ...
                + ncores_data*t_data(t_data > 0)).^2);
        end
        
    end
    
    if flag_metric
        % iteration_number(k) = round(iteration_number(k) / Qc(k)); % only report average iteration number over all sub-problems
        aruntime(k) = aruntime(k)/iteration_number(k);
        atime_l21(k) = atime_l21(k)/iteration_number(k);
        atime_nuclear(k) = atime_nuclear(k)/iteration_number(k);
        atime_master(k) = atime_master(k)/iteration_number(k);
        atime_data(k) = atime_data(k)/iteration_number(k);
        %
        vruntime(k) = (sum_runtime_sqr - iteration_number(k)*aruntime(k)^2)/(iteration_number(k) - 1);
        vtime_l21(k) = (sum_l21_sqr - iteration_number(k)*atime_l21(k)^2)/(iteration_number(k) - 1);
        vtime_nuclear(k) = (sum_nuclear_sqr - iteration_number(k)*atime_nuclear(k)^2)/(iteration_number(k) - 1);
        vtime_master(k) = (sum_master_sqr - iteration_number(k)*atime_master(k)^2)/(iteration_number(k) - 1);
        vtime_data(k) = (sum_data_sqr - iteration_number(k)*atime_data(k)^2)/(iteration_number(k) - 1);  
        % 
        acpu_time(k) = total_cpu_time(k)/iteration_number(k);
        vcpu_time(k) = (sum_cpu_sqr - iteration_number(k)*acpu_time(k)^2)/(iteration_number(k) - 1); 
    
        % compute SNR
        a = SNR(x, x0);
        asnr(k) = mean(a);
        vsnr(k) = var(a);
        a = SNR_log(x, x0);
        asnr_log(k) = mean(a);
        vsnr_log(k) = var(a);
    end
    
    res_ = res_./reshape(operatorNorm, [1, 1, L]); % -> compute operator norm for each channel?
    
    %%
    if flag_display && (Qc(k) == 1) % HyperSARA
        for l = 1:2
            % images
            display_image(x(:,:,chans(l)), clim_log(1,:,l), map_img, fontsize, true);
            export_fig(strcat(savedir,'x_full_hs', num2str(l),'.pdf'), '-transparent','-q101')
            close
            
            % residual images
            display_image(res_(:,:,chans(l)), clim_log(3,:,l), map_img, fontsize, false);
            export_fig(strcat(savedir,'res_full_hs', num2str(l),'.pdf'), '-transparent','-q101')
            close
        end
    end
    
    if flag_display && Qc(k) == 10
        for l = 1:2
            % images
            display_image(x(:,:,chans(l)), clim_log(1,:,l), map_img, fontsize, true);
            export_fig(strcat(savedir,'x_sub', num2str(l),'.pdf'), '-transparent','-q101')
            close
            
            % residual images
            display_image(res_(:,:,chans(l)), clim_log(3,:,l), map_img, fontsize, false);
            export_fig(strcat(savedir,'res_sub', num2str(l),'.pdf'), '-transparent','-q101')
            close
        end
    end
end

%% Display results (table)
if flag_metric
    for k = 1:numel(Qc)
       fprintf("Qc = %i, asnr = %2.2f, std_snr = %1.2e, asnr_log = %2.2f, std_snr_log = %1.2e, iteration_number = %i \n", ...
           Qc(k), asnr(k), sqrt(vsnr(k)), asnr_log(k), sqrt(vsnr_log(k)), iteration_number(k)/Qc(k))
       fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
           aruntime(k), sqrt(vruntime(k)), acpu_time(k), sqrt(vcpu_time(k)));
       fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
           total_runtime(k)/3600, total_cpu_time(k)/3600)
    end

%% Saving results
    save('results_spectral.mat', '-v7.3', 'asnr', 'vsnr', 'asnr_log', ...
        'vsnr_log', 'aruntime', 'vruntime', 'acpu_time', 'vcpu_time', ...
        'iteration_number', ...
        'atime_l21', 'atime_nuclear', 'atime_data', 'atime_master', ...
        'vtime_l21', 'vtime_nuclear', 'vtime_data', 'vtime_master', ...
        'total_runtime', 'total_cpu_time');
end
