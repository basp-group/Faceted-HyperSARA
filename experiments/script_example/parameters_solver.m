% List of default solver-specific parameters (reweighting, pdfb, ellispsoid
% prjection, epsilon update scheme).

%% * Solver
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
param_solver.nu2 = squared_operator_norm_precond;
% regularization parameter nuclear norm
if ~strcmp(algo_version, 'sara')
    param_solver.gamma0 = mu_bar;
end
% regularization parameter l21-norm (soft th parameter)
% ! for SARA, take the value given as an input to the solver
param_solver.gamma = mu;
% id of the cube to be reconstructed (if spectral faceting active)
param_solver.cube_id = ind;
param_solver.backup_frequency = 1;
% print log or not
param_solver.verbose = 2;

% * reweighting
% relative variation (reweighting)
param_solver.reweighting_rel_var = 1e-4;
% minimum number of reweighting iterations
param_solver.reweighting_min_iter = 1;
% TODO: to be removed
param_solver.reweighting_alpha = 1;
param_solver.reweighting_alpha_ff = 1;
% TODO --
% maximum number of reweighting iterations reached (weights updated
% n_reweights times)
param_solver.reweighting_max_iter = max(n_reweights, param_solver.reweighting_min_iter + 1);

% * pdfb
% minimum number of iterations
param_solver.pdfb_min_iter = 10;
% maximum number of iterations
param_solver.pdfb_max_iter = 2000;
% relative variation tolerance
param_solver.pdfb_rel_var = 1e-5;
% tolerance to check data constraints are satisfied
param_solver.pdfb_fidelity_tolerance = 1.01;
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
