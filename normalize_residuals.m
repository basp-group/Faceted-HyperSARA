clc; clear all; close all;
format compact
addpath psf_dr

res_files = ["facethyperSARA_res_it5200_reweight5_gamma5e-06_gamma0_0.01_2b_fouRed2_perc15_adpteps0"; ...
    "facethyperSARA_res_it6400_reweight9_gamma5e-06_gamma0_0.01_2b_fouRed2_perc15_adpteps0";
    "facethyperSARA_res_it9700_reweight20_gamma5e-06_gamma0_0.01_2b_fouRed2_perc15_adpteps0";
    "facethyperSARA_res_it11200_reweight25_gamma5e-06_gamma0_0.01_2b_fouRed2_perc15_adpteps0"];
psf_files = @(channel) ['psf_2b_ind1_16=', num2str(channel)];
channels = [1:17,19,21:32];
rw = [5, 9, 20, 25];

%% Compute all scaling factors
scale_factors = zeros(numel(res_files), 1);
for f = 1:numel(channels)
   psf = fitsread(fullfile("psf_dr", strcat(psf_files(channels(f)), ".fits")));
   scale_factors(f) = max(psf(:));
end
clear psf
save('psf_normalization.mat', '-v7.3', 'scale_factors')

%% Renormalize residual images
scale_factors = reshape(scale_factors, [1, 1, numel(channels)]);
for f = 1:numel(res_files)
    res = fitsread(fullfile(['rw',num2str(rw(f))], strcat(res_files(f), ".fits")));
    res = res./ scale_factors;
    fitswrite(res, ['res_rw', num2str(rw(f)), '.fits'])
end

%% Compute std and kurtosis (only for rw 25)

kurtosis_res = squeeze(kurtosis(res, 0, [1,2]));
akurtosis_res = mean(kurtosis_res);
skurtosis_res = std(kurtosis_res);

std_res = squeeze(std(res, 0, [1,2]));
astd_res = mean(std_res);
sstd_res = std(std_res);

save("metric_rw25.mat", '-v7.3', 'kurtosis_res', 'akurtosis_res', 'skurtosis_res', ...
    'std_res', 'astd_res', 'sstd_res')
