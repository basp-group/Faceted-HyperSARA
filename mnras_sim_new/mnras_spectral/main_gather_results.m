clc; clear all; close all
format compact

addpath ../../lib/faceted-wavelet-transform/src

%% ground truth
image_name = 'cygASband_Cube_256_512_100';
superresolution_factor = 2;
Ny = 256;
Nx = 512;
nChannels = 100;
isnr = 40;
simulation_type = 'spectral';

pathgroundtruth = '../../data/cygASband_Cube_256_512_100.fits';
x0 = fitsread(pathgroundtruth);
N = [size(x0, 1), size(x0, 2)];
nChannels = size(x0, 3);
load('results/cygASband_Cube_256_512_100_spectral/Anorm_hs_Ny=256_Nx=512_L=100.mat', 'operator_norm')
load(strcat('results/cygASband_Cube_256_512_100_spectral/y_', ...
simulation_type,'_',image_name, '_srf=', num2str(superresolution_factor), ...
'_Ny=',num2str(Ny),'_Nx=',num2str(Nx),'_L=', num2str(nChannels), ...
'_snr=', num2str(isnr), ...
'.mat'), 'sigma_noise')

% common parameters
upsilon0 = 3*sigma_noise./sqrt(operator_norm); 

norm2D = @(x) squeeze(sqrt(sum(sum(x.^2, 2), 1)));
SNR_log = @(x, x0, upsilon) 20*log10(norm2D(log10(1 + x0./reshape(upsilon, [1, 1, numel(upsilon)])))./norm2D(log10(1 + x./reshape(upsilon, [1, 1, numel(upsilon)]))-log10(1 + x0./reshape(upsilon, [1, 1, numel(upsilon)]))));

%% SARA
filename = @(ind) fullfile('results/cygASband_Cube_256_512_100_spectral/sara', ...
strcat('spectral_cygASband_Cube_256_512_100_sara_none_srf=2_Ny=256_Nx=512_L=100_Qy=1_Qx=1_Qc=100_ind=', num2str(ind),'_g=3_gb=1_overlap=0_0_hom=0_rwt=heuristic_updreg=0_regtype=heuristic_snr=40_rw=5.mat'));
Q = 1;
ncores_data = 3; % 1 for the data, 2 for the master
ncores_prior = 9;
algo = 'sara';

[x, res, asnr, ssnr, asnr_log, ssnr_log, acpu, scpu, arun, srun, total_cpu_time, total_runtime, iteration_number] = ...
    aggregate_results_spectral(algo, filename, ncores_data, ncores_prior, x0, operator_norm, Q, upsilon0);

a10 = SNR_log(x, x0, 10*upsilon0);
a01 = SNR_log(x, x0, upsilon0/10);

fprintf("snr_log, 10 upsilon: %1.2e (%1.2e) , upsilon/10: %1.2e (%1.2e) \n", mean(a10), std(a10), mean(a01), std(a01));

save(['results_' algo, '_', simulation_type, '.mat'], '-v7.3', 'asnr', 'ssnr', ...
    'asnr_log', 'ssnr_log', 'arun', 'srun', 'acpu', 'scpu', ...
    'iteration_number', 'total_runtime', 'total_cpu_time');

fitswrite(x, "x_sara.fits")
fitswrite(res, "res_sara.fits")

%% HS
filename = @(ind, Qc, alph, alph_bar, ovl) fullfile('results/cygASband_Cube_256_512_100_spectral/hypersara/', ...
['spectral_cygASband_Cube_256_512_100_hypersara_none_srf=2_Ny=256_Nx=512_L=100', ...
'_Qy=1_Qx=1', ...
'_Qc=', num2str(Qc), '_ind=', num2str(ind),'_g=', num2str(alph), '_gb=', num2str(alph_bar), ...
'_overlap=', num2str(ovl),'_', num2str(ovl), ...
'_hom=0_rwt=heuristic3_updreg=0_regtype=heuristic3_snr=40_rw=4.mat']);
alph = 3;
alph_bar = 3;
Qc = [1, 2, 5, 10];
overlap_fraction = 0;
ncores_data = 10;
ncores_prior = 5;
algo = 'hypersara';

for k = 1:numel(Qc)
    qc = Qc(k);
    fprintf("Qc=%i \n", qc)
    filename_ = @(ind) filename(ind, qc, alph, alph_bar, overlap_fraction); 
    [x, res, asnr, ssnr, asnr_log, ssnr_log, acpu, scpu, arun, srun, total_cpu_time, total_runtime, iteration_number] = ...
        aggregate_results_spectral(algo, filename_, ncores_data, ncores_prior, x0, operator_norm, qc, upsilon0);

    a10 = SNR_log(x, x0, 10*upsilon0);
    a01 = SNR_log(x, x0, upsilon0/10);

    fprintf("snr_log, 10 upsilon: %1.2e (%1.2e) , upsilon/10: %1.2e (%1.2e) \n", mean(a10), std(a10), mean(a01), std(a01));

    save(['results_' algo, '_', simulation_type, '.mat'], '-v7.3', 'asnr', 'ssnr', ...
        'asnr_log', 'ssnr_log', 'arun', 'srun', 'acpu', 'scpu', ...
        'iteration_number', 'total_runtime', 'total_cpu_time');
    fitswrite(x, ['x_fhs_Qc=', num2str(qc), '_rw5.fits'])
    fitswrite(res, ['res_fhs_Qc=', num2str(qc), '_rw5.fits'])
end
