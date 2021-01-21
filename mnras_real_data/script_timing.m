clc; 
format compact
addpath ../lib/utils

filename_pattern = @(ind) fullfile('../data', ...
    ['fhs_cw_rd_triangular_Nx=2560_Ny=1536_L=30_Qx=5_Qy=3_Qc=1_overlap=512_ind=',...
     num2str(ind),'_gam0=0.01_gam=5e-06_rw=20.mat']);

files = [2];
nfiles = numel(files);
results_filename = cell(nfiles, 1);
output_filename = 'timing_cyga.mat';

for id = 1:nfiles % 1:nfile
    results_filename{id} = filename_pattern(files(id));
end

[aruntime, vruntime, acpu_time, vcpu_time] = get_timing_facetedHypersara(results_filename, output_filename);


a_ = Qx(k)^2*t_facet(t_facet > 0) + ncores_data*t_data(t_data > 0);
total_cpu_time(k) = sum(a_);
acpu_time(k) = Qx(k)^2*atime_facet(k) + ncores_data*atime_data(k);
vcpu_time(k) = var(a_);