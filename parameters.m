%% Problem size
% TODO: to be filled in

%% NUFFT (gridding parameters)
% oversampling factors for nufft
param_nufft.ox = 2;
% oversampling factors for nufft
param_nufft.oy = 2;
% number of neighbours for nufft
param_nufft.Kx = 8;
% number of neighbours for nufft
param_nufft.Ky = 8;

%% Preconditioning
% number of pixels in the image
param_precond.N = N;
% number of Fourier points (oversampled plane)
param_precond.Nox = ox * Nx;
param_precond.Noy = oy * Ny;
% set weighting type
param_precond.gen_uniform_weight_matrix = 1;
param_precond.uniform_weight_sub_pixels = 1;

%% Blocking
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
% sparam.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi];
% sparam.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi];
% sparam.fpartition = [pi];
% sparam.fpartition = [0, pi];
% sparam.fpartition = [-0.25*pi, 0, 0.25*pi, pi];
% sparam.fpartition = [-64/256*pi, 0, 64/256*pi, pi];

%% * NNLS
% print log or not
param_nnls.verbose = 2;
% stopping criterion
param_nnls.rel_obj = 1e-5;
% max number of iterations
param_nnls.max_iter = 1000;
% saves images at the given iterations
param_nnls.sol_steps = [inf];
param_nnls.beta = 1;

%% Prior

%% Solver
% * general
% estimate of the noise level in SARA space
param_solver.reweighting_sig = sig;
if ~strcmp(algo_version, 'sara')
    % estimate of the noise level in "SVD" spaces
    param_solver.reweighting_sig_bar = sig_bar;
end
% bound on the norm of the Identity operator
param_solver.nu0 = 1;
% bound on the norm of the operator Psi
param_solver.nu1 = 1;
% upper bound on the norm of the measurement operator
param_solver.nu2 = precond_operator_norm; 
param_solver.operator_norm = operator_norm;
% regularization parameter nuclear norm
param_solver.gamma0 = mu_bar;
 % regularization parameter l21-norm (soft th parameter) 
 %! for SARA, take the value given as an input to the solver
param_solver.gamma = mu;
% id of the cube to be reconstructed (if spectral faceting active)
param_solver.cube_id = ind;
param_solver.backup_frequency = 1;
% print log or not
param_solver.verbose = 2;

% * reweighting
% relative variation (reweighting)
param_solver.reweighting_rel_var = 1e-4;
if flag_homotopy
    param_solver.reweighting_alpha = 20;
    param_solver.reweighting_min_iter = 5; % minimum number of reweighting iterations, weights updated reweighting_min_iter times
    param_solver.reweighting_alpha_ff = (1 / param_solver.reweighting_alpha)^(1 / (param_solver.reweighting_min_iter - 1)); 
    % reach the floor level after min_iter updates of the weights
    % 0.63 -> otherwise need 10 reweights minimum
else
    % minimum number of reweighting iterations
    param_solver.reweighting_min_iter = 1; 
    param_solver.reweighting_alpha = 1;
    param_solver.reweighting_alpha_ff = 1;
end
% maximum number of reweighting iterations reached (weights updated 
% nReweights times)
param_solver.reweighting_max_iter = max(nReweights, param_solver.reweighting_min_iter + 1);

% * pdfb
% minimum number of iterations
param_solver.pdfb_min_iter = 10;
% maximum number of iterations
param_solver.pdfb_max_iter = 2000;
% relative variation tolerance
param_solver.pdfb_rel_var = 1e-5;
% tolerance to check data constraints are satisfied
param_solver.pdfb_fidelity_tolerance = 1.01;
param_solver.update_regularization = update_regularization;
param_solver.alph = gam;
param_solver.alph_bar = gam_bar;
% minimum relative variation tolerance (allows stopping earlier if data
% fidelity constraint not about to be satisfied)
param_solver.pdfb_rel_var_low = 5e-6; 

% * ellipsoid projection (if active preconditioning)
% max. number of iterations for the FB algo that implements the 
% preconditioned projection onto the l2 ball min. number of iterations for 
% the FB algo that implements the preconditioned projection onto the l2 
% ball
param_solver.elipse_proj_max_iter = 20;
param_solver.elipse_proj_min_iter = 1;
% precision of the projection onto the ellipsoid
param_solver.elipse_proj_eps = 1e-8; 

% * epsilon update scheme
% flag to activate adaptive epsilon (no need for simulated data)
param_solver.use_adapt_eps = 0; 
% minimum num of iter before stating adjustment
param_solver.adapt_eps_start = 200; 
% tolerance inside the l2 ball
param_solver.adapt_eps_tol_in = 0.99;
% tolerance outside the l2 ball
param_solver.adapt_eps_tol_out = 1.01;
% min num of iter between consecutive updates
param_solver.adapt_eps_steps = 100; 
% bound on the relative change of the solution
param_solver.adapt_eps_rel_var = 5e-5;
% the weight of the update w.r.t the l2 norm of the residual data 
param_solver.adapt_eps_change_percentage = (sqrt(5) - 1) / 2; 
