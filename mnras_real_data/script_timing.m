clc; 
format compact
addpath ../lib/utils

% filename_pattern = @(ind) fullfile('../data', ...
%     ['fhs_cw_rd_triangular_Nx=2560_Ny=1536_L=30_Qx=5_Qy=3_Qc=1_overlap=512_ind=',...
%      num2str(ind),'_gam0=0.01_gam=5e-06_rw=20.mat']);
 
filename_pattern = @(ind) fullfile('/lustre/home/shared/sc004/mnras_faceted_corrected/final_real_data', ...
    ['fhs_cw_rd_triangular_Nx=2560_Ny=1536_L=30_Qx=5_Qy=3_Qc=1_overlap=512_ind=',...
     num2str(ind),'_gam0=0.01_gam=5e-06_rw=20.mat']);

files = [2:5, 7:16];
nfiles = numel(files);
results_filename = cell(nfiles, 1);
% output_filename = '/lustre/home/shared/sc004/mnras_faceted_corrected/final_real_data/timing_faceted_hypersara.mat';
output_filename = 'timing_faceted_hypersara.mat';

for id = 1:nfiles % 1:nfile
    results_filename{id} = filename_pattern(files(id));
end

[aruntime, vruntime, acpu_time, vcpu_time] = get_timing_facetedHypersara(results_filename, output_filename, 16);
