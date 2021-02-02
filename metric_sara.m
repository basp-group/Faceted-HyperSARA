clc; clear all; close all;

%% Compute metric (kurtosis, astd, time) for SARA (with DR)
res = fitsread('res_sara_ddr.fits');
kurtosis_res = squeeze(kurtosis(res, 0, [1,2]));
akurtosis_res = mean(kurtosis_res);
skurtosis_res = std(kurtosis_res);

std_res = squeeze(std(res, 0, [1,2]));
astd_res = mean(std_res);
sstd_res = std(std_res);

skew_res = squeeze(skewness(res,0,[1,2]));
askew_res = mean(skew_res);
sskew_res = std(skew_res);


save("metric_sara.mat", '-v7.3', 'kurtosis_res', 'akurtosis_res', 'skurtosis_res', ...
    'std_res', 'astd_res', 'sstd_res', 'skew_res', 'askew_res', 'sskew_res')

%% Print results
fprintf("astd_res = %1.2e, sstd_res = %1.2e, akurt = %1.2e, skurt = %1.2e, askew = %1.2e, sskew = %1.2e \n", ...
   astd_res, sstd_res, akurtosis_res, skurtosis_res, askew_res, sskew_res)
