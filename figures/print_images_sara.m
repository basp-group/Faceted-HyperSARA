function print_images_sara(results_path, simulation_type, ncores_l11, rw, flag_display, flag_metric)
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
savedir = ['figs_sara_', simulation_type, '/'];
mkdir(savedir)

load(['ground_truth_', simulation_type, '_faceting.mat']) % .mat file generated with 
                                           % `get_ground_truth.m` -> x0, f,
                                           % operatorNorm
[N(1), N(2), L] = size(x0);
x0 = flipud(x0);             % display image using axis convention from ds9
chans = [1, L];

name_pattern = @(ch) fullfile(results_path, ['mnras_faceted_corrected/final_sara_', simulation_type, '/sara_N=', ...
    num2str(N(1)), '_L=',num2str(L),'_ind=1_ch=', num2str(ch),'_gam=0.0001_rw=', num2str(rw),'.mat']);

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

%=========================================================================%
% Plot parameters
%=========================================================================%
clim = [1e-5, 1;  % image
    1e-4, 1;            % error image
    -2.3e-6, 2.3e-6]; % residual images % e-4 before

clim(:,:,2) = [5e-5, 1; % image
    1e-3, 1;                % error image
    -3e-6, 3e-6];     % residual images

fontsize = 25;
map_img = cubehelix(256);

%% Results (HyperSARA, spectrally faceted HyperSARA)

if flag_metric
    x = zeros(N(1), N(2), L);
    
    asnr = 0;
    asnr_log = 0;
    vsnr = 0;
    vsnr_log = 0;
    
    total_runtime = 0;  % assuming all sub-cubes are running in parallel
    total_cpu_time = 0; %
    aruntime = 0;
    vruntime = 0;
    acpu_time = 0;
    vcpu_time = 0;
    atime_l11 = 0; % average time (per iteration)
    atime_data = 0;
    atime_master = 0;
    vtime_l11 = 0; % variance
    vtime_data = 0;
    vtime_master = 0;
    iteration_number = 0;

    sum_runtime_sqr = 0;
    sum_l11_sqr = 0;
    sum_master_sqr = 0;
    sum_data_sqr = 0;
    sum_cpu_sqr = 0;
    
    for l = 1:L
        fileName = name_pattern(l);
        % ncpus =
        load(fileName, 'xsol', 't_l11', 't_master', 't_data', 'end_iter')
        x(:,:,l) = flipud(xsol);
        
        % compute mean and variance over all the files? just 1 for now
        aruntime = aruntime + sum(end_iter(end_iter > 0)); % average runtime per iteration, over all sub-problems
        sum_runtime_sqr = sum_runtime_sqr + sum(end_iter(end_iter > 0).^2);
        
        atime_l11 = atime_l11 + sum(t_l11(t_l11 > 0));
        atime_master = atime_master + sum(t_master(t_master > 0));
        atime_data = atime_data + sum(t_data(t_data > 0));
        sum_l11_sqr = sum_l11_sqr + sum(t_l11(t_l11 > 0).^2);
        sum_master_sqr = sum_master_sqr + sum(t_master(t_master > 0).^2);
        sum_data_sqr = sum_data_sqr + sum(t_data(t_data > 0).^2);
        
        % average number of iterations over all sub-problems
        iteration_number = iteration_number + sum(end_iter > 0);
        total_runtime = max(total_runtime, sum(end_iter(end_iter > 0)));
        total_cpu_time = total_cpu_time + ncores_l11*sum(t_l11(t_l11 > 0)) + ...
            + sum(t_master(t_master > 0)) ...
            + sum(t_data(t_data > 0));
        sum_cpu_sqr = sum_cpu_sqr + ncores_l11*sum((t_l11(t_l11 > 0) + t_master(t_master > 0) ...
            + t_data(t_data > 0)).^2);
    end
    
    aruntime = aruntime/iteration_number;
    atime_l11 = atime_l11/iteration_number;
    atime_master = atime_master/iteration_number;
    atime_data = atime_data/iteration_number;
    %
    vruntime = (sum_runtime_sqr - iteration_number*aruntime^2)/(iteration_number - 1);
    vtime_l11 = (sum_l11_sqr - iteration_number*atime_l11^2)/(iteration_number - 1);
    vtime_master = (sum_master_sqr - iteration_number*atime_master^2)/(iteration_number - 1);
    vtime_data = (sum_data_sqr - iteration_number*atime_data^2)/(iteration_number - 1);
    %
    acpu_time = total_cpu_time/iteration_number;
    vcpu_time = (sum_cpu_sqr - iteration_number*acpu_time^2)/(iteration_number - 1);
    
    % compute SNR
    a = SNR(x, x0);
    asnr = mean(a);
    vsnr = var(a);
    a = SNR_log(x, x0);
    asnr_log = mean(a);
    vsnr_log = var(a);
    
    %% Display results (table)
    fprintf("SARA %s, asnr = %2.2f, std_snr = %1.2e, asnr_log = %2.2f, std_snr_log = %1.2e, iteration_number = %i \n", ...
       simulation_type, asnr, sqrt(vsnr), asnr_log, sqrt(vsnr_log), iteration_number)
    fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
       aruntime, sqrt(vruntime), acpu_time, sqrt(vcpu_time));
    fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
       total_runtime/3600, total_cpu_time/3600)

    %% Saving results
    save(['results_sara_', simulation_type, '.mat'], '-v7.3', 'asnr', 'vsnr', 'asnr_log', ...
        'vsnr_log', 'aruntime', 'vruntime', 'acpu_time', 'vcpu_time', ...
        'iteration_number', ...
        'atime_l11', 'atime_data', 'atime_master', ...
        'vtime_l11', 'vtime_data', 'vtime_master', ...
        'total_runtime', 'total_cpu_time');
end
    
%%
if flag_display
    for l = 1:2
        fileName = name_pattern(chans(l));
        load(fileName, 'xsol', 'res');

        xsol = flipud(xsol);
        res = flipud(res)/operatorNorm(chans(l));

        % images
        display_image(xsol, clim(1,:,l), map_img, fontsize, true);
        export_fig(strcat(savedir,'x_sara', num2str(l),'.pdf'), '-transparent','-q101')
        close

        % residual images
        display_image(res, clim(3,:,l), map_img, fontsize, false);
        export_fig(strcat(savedir,'res_sara', num2str(l),'.pdf'), '-transparent','-q101')
        close
    end
end
