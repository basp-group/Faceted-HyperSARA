% List of default parameters related to the problem configuration (nufft,
% preconditioning, blocking, NNLS, SARA dictionary)

%% * Problem size
% TODO: to be completed for real data, revise blocking for real data
speed_of_light = 299792458;
% seed for the random number generation (used only when generating data)
seed = 1;

%% * NUFFT (gridding parameters)
% oversampling factor (x)
param_nufft.ox = 2;
% oversampling factor (y)
param_nufft.oy = 2;
% number of neighbour (x)
param_nufft.Kx = 7;
% number of neighbour (y)
param_nufft.Ky = 7;

%% * Preconditioning
% number of pixels in the image
param_precond.N = Nx * Ny;
% number of Fourier points (oversampled plane)
param_precond.Nox = param_nufft.ox * Nx;
param_precond.Noy = param_nufft.oy * Ny;
% set weighting type
param_precond.gen_uniform_weight_matrix = 1;
param_precond.uniform_weight_sub_pixels = 1;

%% * Blocking
% density-based
param_blocking.use_density_partitioning = 0;
param_blocking.density_partitioning_no = 1;
% uniform
param_blocking.use_uniform_partitioning = 0;
param_blocking.uniform_partitioning_no = 4;
% equal-size
param_blocking.use_equal_partitioning = 1;
param_blocking.equal_partitioning_no = 1;
% manual
param_blocking.use_manual_partitioning = 0;
param_blocking.use_manual_frequency_partitioning = 0;
% partition (symetrically) of the data to nodes (frequency ranges)
param_blocking.fpartition = [icdf('norm', 0.25, 0, pi / 4), 0, icdf('norm', 0.75, 0, pi / 4), pi];

%% * Prior
% depth of the wavelet decompositions
nlevel = 4;
% wavelet dictionaries
% ! always specify Dirac basis ('self') in last position if ever used
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};
% length of the filters (0 corresponding to the 'self' basis)
filter_length = [2 * (1:8)'; 0];
