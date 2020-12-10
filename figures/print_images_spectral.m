function print_images_spectral(results_path)
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
% doComputeMetric = true;                      % save metrics in a .txt file
display = true;

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
col = @(x) x(:);
norm2D = @(x) squeeze(sqrt(sum(sum(x.^2, 2), 1)));
norm2 = @(x) squeeze(sqrt(sum(sum(x(:).^2, 2), 1)));
DR = @(x, res, n) sqrt(prod(N))*squeeze(max(max(x,[],2),[],1))*n./norm2D(res); % per band input
SNR = @(x, x0) 20*log10(norm2D(x0)./norm2D(x - x0));
% SNR_log = @(x, x0) SNR(log10(x+eps), log10(x0+eps));
SNR_log = @(x, x0) 20*log10(norm2D(log10(x0+eps))./norm2D(log10(x+eps)-log10(x0+eps)));

asnr = zeros(numel(Qc), 1);
asnr_log = zeros(numel(Qc), 1);
vsnr = zeros(numel(Qc), 1);
vsnr_log = zeros(numel(Qc), 1);
runtime = zeros(numel(Qc), 1);
cpu_time = zeros(numel(Qc), 1);
atime_l21 = zeros(numel(Qc), 1); % average time (per iteration)
atime_nuclear = zeros(numel(Qc), 1);
atime_data = zeros(numel(Qc), 1);
atime_master = zeros(numel(Qc), 1);
vtime_l21 = zeros(numel(Qc), 1); % variance
vtime_nuclear = zeros(numel(Qc), 1);
vtime_data = zeros(numel(Qc), 1);
vtime_master = zeros(numel(Qc), 1);

%=========================================================================%
% Plot parameters
%=========================================================================%
clim_log = [-5, 0;  % image
    -4, 0;            % error image
    -3.5e-6, 3.5e-6]; % residual images % e-4 before

clim_log(:,:,2) = [-4, 0; % image
    -3, 0;                % error image
    -3.5e-6, 3.5e-6];     % residual images

fontsize = 25;
map_img = cubehelix(256);

%% Display ground truth image (spectral faceting)
for l = 1:2
    % ground truth
    display_image(log10(x0(:,:,chans(l))), clim_log(1,:,l), map_img, fontsize);
    export_fig(strcat(savedir,'x', num2str(l),'_sub.pdf'), '-transparent','-q101')
    close
end

%% Results (HyperSARA, spectrally faceted HyperSARA)
for k = 1:numel(Qc)
    x = zeros(N(1), N(2), L);
    res_ = zeros(N(1), N(2), L);
    channels = split_range_interleaved(Qc(k), L);
    for ind = 1:Qc(k)
        fileName = name_pattern(Qc(k), ind);
        % ncpus =
        load(fileName, 'xsol', 'res', 't_l21', 't_nuclear', 't_master', 't_data', 'end_iter') % , 'param'
        % operator_norm = param.nu2;
        x(:,:,channels{ind}) = flipud(xsol);
        res_(:,:,channels{ind}) = flipud(res); 
        
        % compute mean and variance over all the files? just 1 for now
        if ind == 1
            runtime(k) = mean(end_iter(end_iter > 0)); % average runtime per iteration
            atime_l21(k) = mean(t_l21(t_l21 > 0));
            atime_nuclear(k) = mean(t_nuclear(t_nuclear > 0));
            atime_master(k) = mean(t_master(t_master > 0));
            atime_data(k) = mean(t_data(t_data > 0));
            vtime_l21(k) = var(t_l21(t_l21 > 0));
            vtime_nuclear(k) = var(t_nuclear(t_nuclear > 0));
            vtime_master(k) = var(t_master(t_master > 0));
            vtime_data(k) = var(t_data(t_data > 0));
            % average cpu time per iteration
            cpu_time(k) = atime_l21(k) + atime_nuclear(k) + ...
                atime_master(k) + numel(channels{ind})*atime_data(k);  
        end
    end
    
    % compute SNR
    a = SNR(x, x0);
    asnr(k) = mean(a);
    vsnr(k) = var(a);
    a = SNR_log(x, x0);
    asnr_log(k) = mean(a);
    vsnr_log(k) = var(a);
    
    res_ = res_./reshape(operatorNorm, [1, 1, L]); % -> compute operator norm for each channel?
    
    %%
    if display && (Qc(k) == 1) % HyperSARA
        for l = 1:2
            % images
            display_image(log10(x(:,:,chans(l))), clim_log(1,:,l), map_img, fontsize);
            export_fig(strcat(savedir,'x_full_hs', num2str(l),'.pdf'), '-transparent','-q101')
            close
            
            % residual images
            display_image(res_(:,:,chans(l)), clim_log(3,:,l), map_img, fontsize);
            export_fig(strcat(savedir,'res_full_hs', num2str(l),'.pdf'), '-transparent','-q101')
            close
        end
    end
    
    if display && Qc(k) == 10
        for l = 1:2
            % images
            display_image(log10(x(:,:,chans(l))), clim_log(1,:,l), map_img, fontsize);
            export_fig(strcat(savedir,'x_sub', num2str(l),'.pdf'), '-transparent','-q101')
            close
            
            % residual images
            display_image(res_(:,:,chans(l)), clim_log(3,:,l), map_img, fontsize);
            export_fig(strcat(savedir,'res_sub', num2str(l),'.pdf'), '-transparent','-q101')
            close
        end
    end
end

%% Display results (table)
for k = 1:numel(Qc)
   fprintf("Qc = %i, asnr = %2.2f, vsnr = %1.2e, asnr_log = %2.2f, vsnr_log = %1.2e, runtime = %.2f, cpu_time = %.2f \n", ...
       Qc(k), asnr(k), vsnr(k), asnr_log(k), vsnr_log(k), runtime(k), cpu_time(k))
end

%% Saving results
save('results_spectral.mat', '-v7.3', 'asnr', 'vsnr', 'asnr_log', ...
    'vsnr_log', 'runtime', 'cpu_time', ...
    'atime_l21', 'atime_nuclear', 'atime_data', 'atime_master', ...
    'vtime_l21', 'vtime_nuclear', 'vtime_data', 'vtime_master')
