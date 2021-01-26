clc; 
format compact
addpath ../lib/utils

% filename_pattern = @(ind) fullfile('../data', ...
%     ['fhs_cw_rd_triangular_Nx=2560_Ny=1536_L=30_Qx=5_Qy=3_Qc=1_overlap=512_ind=',...
%      num2str(ind),'_gam0=0.01_gam=5e-06_rw=20.mat']);

%% Faceted HyperSARA
res_path = '/lustre/home/shared/sc004/PURE_MJ/new_results/final_real_data/facethyperSARA_res_norm_it6300_reweight20_gamma5e-06_gamma0_0.01_2b_fouRed2_perc15_adpteps0.fits';
results_filename = "/lustre/home/shared/sc004/adrianj/new_results/facethyperSARA_dr_co_w_real_1_16_5e-06_20.mat";
nfacets = 15; % Qy = 3, Qx = 5, ncores_data = 15
output_filename = 'timing_faceted_hypersara.mat';

[aruntime, vruntime, acpu_time, vcpu_time, total_runtime, total_cpu_time, iteration_number] = get_timing_facetedHypersara(results_filename, nfacets);

fprintf("HyperSARA: Qc = %i, iteration_number = %i \n", ...
        Qc, iteration_number/Qc)
fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
    aruntime, sqrt(vruntime), acpu_time, sqrt(vcpu_time));
fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
    total_runtime/3600, total_cpu_time/3600)

% load full (normalised) residual cube
res = fitsread(res_path);

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

save("metric_fhs_dr.mat", '-v7.3', 'kurtosis_res', 'akurtosis_res', 'skurtosis_res', ...
    'std_res', 'astd_res', 'sstd_res', 'skew_res', 'askew_res', 'sskew_res', ...
    'aruntime', 'vruntime', 'acpu_time', 'vcpu_time', ...
    'iteration_number', ...
    'total_runtime', 'total_cpu_time')

% Print results
fprintf("astd_res = %1.2e, sstd_res = %1.2e, akurt = %1.2e, skurt = %1.2e, askew = %1.2e, sskew = %1.2e \n", ...
   astd_res, sstd_res, akurtosis_res, skurtosis_res, askew_res, sskew_res)

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
fprintf("SARA: astd_res = %1.2e, sstd_res = %1.2e, akurt = %1.2e, skurt = %1.2e, askew = %1.2e, sskew = %1.2e \n", ...
   astd_res, sstd_res, akurtosis_res, skurtosis_res, askew_res, sskew_res)

save("metric_sara_dr.mat", '-v7.3', 'kurtosis_res', 'akurtosis_res', 'skurtosis_res', ...
   'std_res', 'astd_res', 'sstd_res', 'skew_res', 'askew_res', 'sskew_res')
%    , ...
%    'aruntime', 'vruntime', 'acpu_time', 'vcpu_time', ...
%    'iteration_number', ...
%    'total_runtime', 'total_cpu_time')

%! to be modified
% filename_pattern = @(ind) fullfile('/lustre/home/shared/sc004/mnras_faceted_corrected/final_real_data', ...
%     ['fhs_cw_rd_triangular_Nx=2560_Ny=1536_L=30_Qx=5_Qy=3_Qc=1_overlap=512_ind=',...
%      num2str(ind),'_gam0=0.01_gam=5e-06_rw=20.mat']);
% for id = 1:Qc % 1:nfile
%     results_filename{id} = filename_pattern(files(id));
% end
% output_filename = 'timing_faceted_sara.mat';
%
% [aruntime, vruntime, acpu_time, vcpu_time, total_runtime, total_cpu_time, iteration_number] = get_timing_sara(results_filename);
% 
% if L < ... % if all sub-cubes not available
%     % estimate total cpu time based on number of iterations and average
%     % time per iteration per subcube
%     tcpu = total_cpu_time + (...-L)*(iteration_number/Qc)*acpu_time;
%       
%     fprintf("Qc = %i, iteration_number = %i \n", ...
%         Qc, iteration_number/Qc)
%     fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
%         aruntime, sqrt(vruntime), acpu_time, sqrt(vcpu_time));
%     fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
%         total_runtime/3600, tcpu/3600)
% end
