clc; 
format compact
addpath ../lib/utils

% filename_pattern = @(ind) fullfile('../data', ...
%     ['fhs_cw_rd_triangular_Nx=2560_Ny=1536_L=30_Qx=5_Qy=3_Qc=1_overlap=512_ind=',...
%      num2str(ind),'_gam0=0.01_gam=5e-06_rw=20.mat']);

%% Faceted HyperSARA
filename_pattern = @(ind) fullfile('/lustre/home/shared/sc004/mnras_faceted_corrected/final_real_data', ...
    ['fhs_cw_rd_triangular_Nx=2560_Ny=1536_L=30_Qx=5_Qy=3_Qc=1_overlap=512_ind=',...
     num2str(ind),'_gam0=0.01_gam=5e-06_rw=20.mat']);

nfacets = 15; % Qy = 3, Qx = 5, ncores_data = 15
files = [2:5, 7:16]; % [1:16];
Qc = numel(files);   % number of spectral sub-cubes
results_filename = cell(Qc, 1);
% output_filename = '/lustre/home/shared/sc004/mnras_faceted_corrected/final_real_data/timing_faceted_hypersara.mat';
output_filename = 'timing_faceted_hypersara.mat';

for id = 1:Qc % 1:nfile
    results_filename{id} = filename_pattern(files(id));
end

[aruntime, vruntime, acpu_time, vcpu_time, total_runtime, total_cpu_time, ... 
    iteration_number] = get_timing_facetedHypersara(results_filename, output_filename, nfacets);

if Qc < 16 % if all sub-cubes not available
    % estimate total cpu time based on number of iterations and average
    % time per iteration per subcube
    tcpu = total_cpu_time + (16-Qc)*(iteration_number/Qc)*acpu_time;
      
    fprintf("Qc = %i, iteration_number = %i \n", ...
        Qc, iteration_number/Qc)
    fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
        aruntime, sqrt(vruntime), acpu_time, sqrt(vcpu_time));
    fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
        total_runtime/3600, tcpu/3600)
end

%% SARA (only placeholders for now)

%! to be modified
% filename_pattern = @(ind) fullfile('/lustre/home/shared/sc004/mnras_faceted_corrected/final_real_data', ...
%     ['fhs_cw_rd_triangular_Nx=2560_Ny=1536_L=30_Qx=5_Qy=3_Qc=1_overlap=512_ind=',...
%      num2str(ind),'_gam0=0.01_gam=5e-06_rw=20.mat']);
% for id = 1:Qc % 1:nfile
%     results_filename{id} = filename_pattern(files(id));
% end
% output_filename = 'timing_faceted_sara.mat';
%
% [aruntime, vruntime, acpu_time, vcpu_time, total_runtime, total_cpu_time, iteration_number] = get_timing_sara(results_filename, metric_filename);
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
