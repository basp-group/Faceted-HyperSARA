% List of default parameters related to the problem configuration (nufft,
% preconditioning, blocking, NNLS, SARA dictionary)

%% * Problem size
% TODO: to be completed

%% * NUFFT (gridding parameters)
param_nufft.ox = 2; % oversampling factor (x)
param_nufft.oy = 2; % oversampling factor (y)
param_nufft.Kx = 7; % 8% number of neighbour (x)
param_nufft.Ky = 7; % 8% number of neighbour (y)

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

%% * NNLS (estimation of the l2 contraint)
% flag to activate NNLS
generate_eps_nnls = false;
% original image, used to compute the SNR
% param_nnls.im = im;
% print log or not
param_nnls.verbose = 2;
% stopping criterion
param_nnls.rel_obj = 5e-4; % 1e-5 is too low !!
% max number of iterations
param_nnls.max_iter = 1000;
% saves images at the given iterations
param_nnls.sol_steps = [inf];
param_nnls.beta = 1;
