clc; 
format compact
addpath ../../lib/utils


%% Faceted HyperSARA
% res_path = '/lustre/home/shared/sc004/PURE_MJ/new_results/final_real_data/facethyperSARA_res_norm_it6300_reweight20_gamma5e-06_gamma0_0.01_2b_fouRed2_perc15_adpteps0.fits';
% results_filename = "/lustre/home/shared/sc004/adrianj/new_results/facethyperSARA_dr_co_w_real_1_16_5e-06_20.mat";
res_path = '../mnras_faceted_corrected/final_dr/r_fhs_dr.fits';
results_filename = "../mnras_faceted_corrected/final_dr/gamma5e-6/facethyperSARA_dr_co_w_real_1_16_5e-06_0.001_15_adpteps0_22.mat";        
            
nfacets = 15; % Qy = 3, Qx = 5, ncores_data = 15
ncores_data = 15;
output_filename = 'timing_faceted_hypersara.mat';

[aruntime, vruntime, acpu_time, vcpu_time, total_runtime, total_cpu_time, iteration_number] = get_timing_facetedHypersara(results_filename, nfacets, ncores_data);

fprintf("HyperSARA: iteration_number = %i \n", ...
        iteration_number)
fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
    aruntime, sqrt(vruntime), acpu_time, sqrt(vcpu_time));
fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
    total_runtime/3600, total_cpu_time/3600)

% load full (normalised) residual cube
res = fitsread(res_path);
ares = mean(res, 3);
n_channels = size(res, 3);

std_ares = std(ares(:));

kurtosis_res = squeeze(kurtosis(res, 0, [1,2]));
akurtosis_res = mean(kurtosis_res);
skurtosis_res = std(kurtosis_res);

std_res = squeeze(std(res, 0, [1,2]));
astd_res = mean(std_res);
sstd_res = std(std_res);

skew_res = squeeze(skewness(res,0,[1,2]));
askew_res = mean(skew_res);
sskew_res = std(skew_res);
clear res ares

save("metric_fhs_dr.mat", '-v7.3', 'kurtosis_res', 'akurtosis_res', 'skurtosis_res', ...
    'std_res', 'astd_res', 'sstd_res', 'std_ares', 'skew_res', 'askew_res', 'sskew_res', ...
    'aruntime', 'vruntime', 'acpu_time', 'vcpu_time', ...
    'iteration_number', ...
    'total_runtime', 'total_cpu_time')

% Print results
fprintf("std_res [1st, last] = [%1.2e, %1.2e], astd_res = %1.2e, sstd_res = %1.2e, std_ares = %1.2e, akurt = %1.2e, skurt = %1.2e, askew = %1.2e, sskew = %1.2e \n", ...
   std_res(1), std_res(end), astd_res, sstd_res, std_ares, akurtosis_res, skurtosis_res, askew_res, sskew_res)

% This part was only needed when interpolating results (missing cubes)
% if Qc < 16 % if all sub-cubes not available
%     % estimate total cpu time based on number of iterations and average
%     % time per iteration per subcube
%     tcpu = total_cpu_time + (16-Qc)*(iteration_number/Qc)*acpu_time;
      
%     fprintf("Qc = %i, iteration_number = %i \n", ...
%         Qc, iteration_number/Qc)
%     fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
%         aruntime, sqrt(vruntime), acpu_time, sqrt(vcpu_time));
%     fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
%         total_runtime/3600, tcpu/3600)
% end

%% SARA (only placeholders for now)

% load full (normalised) residual cube
res_path = '/lustre/home/shared/sc004/PURE_MJ/res_sara_ddr.fits';
res = fitsread(res_path);
ares = mean(res, 3);

std_ares = std(ares(:));

kurtosis_res = squeeze(kurtosis(res, 0, [1,2]));
akurtosis_res = mean(kurtosis_res);
skurtosis_res = std(kurtosis_res);

std_res = squeeze(std(res, 0, [1,2]));
astd_res = mean(std_res);
sstd_res = std(std_res);

skew_res = squeeze(skewness(res,0,[1,2]));
askew_res = mean(skew_res);
sskew_res = std(skew_res);
clear res

% Print results
fprintf("std_res [1st, last] = [%1.2e, %1.2e], astd_res = %1.2e, sstd_res = %1.2e, std_ares = %1.2e, akurt = %1.2e, skurt = %1.2e, askew = %1.2e, sskew = %1.2e \n", ...
   std_res(1), std_res(end), astd_res, sstd_res, std_ares, akurtosis_res, skurtosis_res, askew_res, sskew_res)

%% ! to be modified
n_available_channels = 1;
n_channels = 30; % total number of channels in the cube
% filename_pattern = @(ind) fullfile('/lustre/home/shared/sc004/PURE_MJ/monoSARA', ...
%     ['SARA_ch', num2str(ind), '_prec_reweight20_rw_alpha0.028147_gamma5e-06_2b_fouRed2_perc15']);

filename_pattern = @(ind) fullfile('../mnras_faceted_corrected/final_dr/monoSARA', ...
    ['SARA_ch', num2str(ind), '_prec_reweight20_rw_alpha0.028147_gamma5e-06_2b_fouRed2_perc15.mat']);

results_filename = strings(n_available_channels, 1);
for id = 1:n_available_channels
    results_filename(id) = filename_pattern(id);
end
output_filename = 'timing_faceted_sara.mat';

%%
[aruntime, vruntime, acpu_time, vcpu_time, total_runtime, total_cpu_time, iteration_number] = get_timing_sara(results_filename);

if n_available_channels < n_channels % if all channels not available
    % estimate total cpu time based on number of iterations and average
    % time per iteration per subcube
    tcpu = total_cpu_time + (n_channels-n_available_channels)*(iteration_number/n_available_channels)*acpu_time;
     
    fprintf("SARA, iteration_number = %i \n", ...
        iteration_number)
    fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
        aruntime, sqrt(vruntime), acpu_time, sqrt(vcpu_time));
    fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
        total_runtime/3600, tcpu/3600)
end


% Print results
fprintf("SARA: std_res [1st, last] = [%1.2e, %1.2e], astd_res = %1.2e, sstd_res = %1.2e, std_ares = %1.2e, akurt = %1.2e, skurt = %1.2e, askew = %1.2e, sskew = %1.2e \n", ...
   std_res(1), std_res(end), std_astd_res, sstd_res, strd_ares, akurtosis_res, skurtosis_res, askew_res, sskew_res)

%%
% taking into account the reference 10h from the previous submission (consistency reasons)
iteration_number = 10000; % from log file of previous run
tcpu = acpu_time*iteration_number + (n_channels-n_available_channels)*(iteration_number/n_available_channels)*acpu_time;     
fprintf("SARA, iteration_number = %i \n", ...
    iteration_number)
fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
    aruntime, sqrt(vruntime), acpu_time, sqrt(vcpu_time));
fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
    aruntime*iteration_number/3600, tcpu/3600)

save("metric_sara_dr.mat", '-v7.3', 'kurtosis_res', 'akurtosis_res', 'skurtosis_res', ...
'std_res', 'astd_res', 'sstd_res', 'std_ares', 'skew_res', 'askew_res', 'sskew_res', ...
'aruntime', 'vruntime', 'acpu_time', 'vcpu_time', ...
'iteration_number', ...
'total_runtime', 'total_cpu_time')

%% Extract timing from log file (manually)
tab = load("timing_dr.txt"); % iter, runtime, t_facet, t_data, rel_var

nfacets = 15; % Qy = 3, Qx = 5, ncores_data = 15
ncores_data = 15;

total_runtime = sum(tab(:,2));
aruntime = mean(tab(:,2));
std_runtime = std(tab(:,2));
atime_facet = mean(tab(:,3));
std_facet = std(tab(:,3));
atime_data = mean(tab(:,4));
std_data = std(tab(:,4));

% average number of iterations over all sub-problems
iteration_number = size(tab,1);
cpu_time = nfacets*tab(:,3) + ncores_data*tab(:,4);

total_cpu_time = nfacets*sum(tab(:,3)) ...
    + ncores_data*sum(tab(:,4));
% sum_cpu_sqr = sum_cpu_sqr + sum((ncores_facets*tab(:,3) ...
%     + ncores_data*tab(:,4)).^2);
acpu_time = mean(cpu_time);
std_cpu_time = std(cpu_time);

fprintf("iteration_number = %i \n", ...
    iteration_number)
fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
    aruntime, std_runtime, acpu_time, std_cpu_time);
fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
    total_runtime/3600, total_cpu_time/3600)
