clc; clear all; close all
format compact

%% ground truth
image_name = 'cygASband_Cube_1024_2048_20';
superresolution_factor = 2;
Ny = 1024;
Nx = 2048;
nChannels = 20;
isnr = 40;
simulation_type = 'spatial';

pathgroundtruth = '../../data/cygASband_Cube_1024_2048_20.fits';
x0 = fitsread(pathgroundtruth);
N = [size(x0, 1), size(x0, 2)];
nChannels = size(x0, 3);
load('results/cygASband_Cube_1024_2048_20_spatial/Anorm_hs_Ny=1024_Nx=2048_L=20.mat', 'operator_norm')
load(strcat('results/cygASband_Cube_1024_2048_20_spatial/y_', ...
simulation_type,'_',image_name, '_srf=', num2str(superresolution_factor), ...
'_Ny=',num2str(Ny),'_Nx=',num2str(Nx),'_L=', num2str(nChannels), ...
'_snr=', num2str(isnr), ...
'.mat'), 'sigma_noise')

% common parameters
upsilon0 = 3*sigma_noise./sqrt(operator_norm); 

norm2D = @(x) squeeze(sqrt(sum(sum(x.^2, 2), 1)));
SNR_log = @(x, x0, upsilon) 20*log10(norm2D(log10(1 + x0./reshape(upsilon, [1, 1, numel(upsilon)])))./norm2D(log10(1 + x./reshape(upsilon, [1, 1, numel(upsilon)]))-log10(1 + x0./reshape(upsilon, [1, 1, numel(upsilon)]))));

%% SARA
filename = @(ind) fullfile('results/cygASband_Cube_1024_2048_20_spatial/sara', ...
strcat('spatial_cygASband_Cube_1024_2048_20_sara_none_srf=2_Ny=1024_Nx=2048_L=20_Qy=1_Qx=1_Qc=20_ind=', num2str(ind),'_g=3_gb=1_overlap=0_0_hom=0_rwt=heuristic_updreg=0_regtype=heuristic_snr=40_rw=5.mat'));
Q = 1;
ncores_data = 3; % 1 for the data, 2 for the master
ncores_prior = 9;
algo = 'sara';

[x, res, asnr, ssnr, asnr_log, ssnr_log, acpu, scpu, arun, srun, total_cpu_time, total_runtime, iteration_number] = ...
    aggregate_results(algo, filename, ncores_data, ncores_prior, x0, operator_norm, Q, upsilon0);

a10 = SNR_log(x, x0, 10*upsilon0);
a01 = SNR_log(x, x0, upsilon0/10);

fprintf("snr_log, 10 upsilon: %1.2e (%1.2e) , upsilon/10: %1.2e (%1.2e) \n", mean(a10), std(a10), mean(a01), std(a01));

save(['results_' algo, '_', simulation_type, '.mat'], '-v7.3', 'asnr', 'ssnr', ...
    'asnr_log', 'ssnr_log', 'arun', 'srun', 'acpu', 'scpu', ...
    'iteration_number', 'total_runtime', 'total_cpu_time');

fitswrite(x, "x_sara.fits")
fitswrite(res, "res_sara.fits")


%% HS
filename = @(Q, alph, alph_bar, ovl) fullfile('results/cygASband_Cube_1024_2048_20_spatial/hypersara/heuristic', ...
['spatial_cygASband_Cube_1024_2048_20_hypersara_none_srf=2_Ny=1024_Nx=2048_L=20', ...
'_Qy=', num2str(Q), '_Qx=', num2str(Q), ...
'_Qc=1_ind=1_g=', num2str(alph), '_gb=', num2str(alph_bar), ...
'_overlap=', num2str(ovl),'_', num2str(ovl), ...
'_hom=0_rwt=heuristic2_updreg=0_regtype=heuristic2_snr=40_rw=5.mat']);
alph = 1;
alph_bar = 3;
Q = 1;
overlap_fraction = 0;
ncores_data = 20;
ncores_prior = 2;
algo = 'hypersara';

[x, res, asnr, ssnr, asnr_log, ssnr_log, acpu, scpu, arun, srun, total_cpu_time, total_runtime, iteration_number] = ...
    aggregate_results(algo, filename, ncores_data, ncores_prior, x0, operator_norm, Q, upsilon0);

a10 = SNR_log(x, x0, 10*upsilon0);
a01 = SNR_log(x, x0, upsilon0/10);

fprintf("snr_log, 10 upsilon: %1.2e (%1.2e) , upsilon/10: %1.2e (%1.2e) \n", mean(a10), std(a10), mean(a01), std(a01));

save(['results_' algo, '_', simulation_type, '.mat'], '-v7.3', 'asnr', 'ssnr', ...
    'asnr_log', 'ssnr_log', 'arun', 'srun', 'acpu', 'scpu', ...
    'iteration_number', 'total_runtime', 'total_cpu_time');
fitswrite(x, "x_hs.fits")
fitswrite(res, "res_hs.fits")


%% FHS
filename = @(Q, alph, alph_bar, ovl) fullfile('results/cygASband_Cube_1024_2048_20_spatial/cw/', ...
['spatial_cygASband_Cube_1024_2048_20_cw_triangular_srf=2_Ny=1024_Nx=2048_L=20', ...
'_Qy=', num2str(Q), '_Qx=', num2str(Q), ...
'_Qc=1_ind=1_g=', num2str(alph), '_gb=', num2str(alph_bar), ...
'_overlap=', num2str(ovl),'_', num2str(ovl), ...
'_hom=0_rwt=heuristic2_updreg=0_regtype=heuristic2_snr=40_rw=5.mat']);
alph = 1;
alph_bar = 3;
Q = [4];
overlap_fraction = [0, 0.1, 0.25, 0.4];
ncores_data = 20;
ncores_prior = Q.^2;
algo = 'cw';

for k = 1:numel(Q)
    q = Q(k);
    for l = 1:numel(overlap_fraction)
        ovl = overlap_fraction(l);
        fprintf("Q=%i, ovl=%1.2e \n", q, ovl)
        [x, res, asnr, ssnr, asnr_log, ssnr_log, acpu, scpu, arun, srun, total_cpu_time, total_runtime, iteration_number] = ...
            aggregate_results(algo, filename(q, alph, alph_bar, ovl), ncores_data, ncores_prior, x0, operator_norm, q, upsilon0);

        a10 = SNR_log(x, x0, 10*upsilon0);
        a01 = SNR_log(x, x0, upsilon0/10);

        fprintf("snr_log, 10 upsilon: %1.2e (%1.2e) , upsilon/10: %1.2e (%1.2e) \n", mean(a10), std(a10), mean(a01), std(a01));

        save(['results_' algo, '_', simulation_type, '.mat'], '-v7.3', 'asnr', 'ssnr', ...
            'asnr_log', 'ssnr_log', 'arun', 'srun', 'acpu', 'scpu', ...
            'iteration_number', 'total_runtime', 'total_cpu_time');
        fitswrite(x, ['x_fhs_Q=', num2str(q), '_ovl=', num2str(ovl), '.fits'])
        fitswrite(res, ['res_fhs_Q=', num2str(q), '_ovl=', num2str(ovl), '.fits'])
    end
end
