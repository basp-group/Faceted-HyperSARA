function [param_global, param_solver, param_nufft, ...
   param_precond, param_wproj, dict] = ...
    read_json_configuration(json_filename)
% Read algorithmic configuration parameters defined in an input ``.json``
% file.
%
% Parameters
% ----------
% json_filename : string
%     Name of the .json configuration file.
% Returns
% -------
% param_global : struct
%     Structure containing global configuration parameters.
% param_solver : struct
%     Structure containing solver parameters (reweighting scheme, PDFB,
%     projection onto the ellipsoid, noise estimation on the fly).
% param_nufft : struct
%     Structure containing the non-uniform FFT parameters.
% param_precond : struct
%     Structure containing parameters used to define the preconditioning
%     matrix involved in the PDFB algorithm.
% param_wproj : struct
%     Structure containing `w`-projection parameters.
% dict : struct
%     Structure defining the SARA dictionary.
%
% Important
% ---------
% The ``.json`` file should contain fields corresponding to the following structure attributes.
%
% param_solver.reweighting_rel_var  (double)
%     Relative variation stopping criterion (reweighting).
% param_solver.reweighting_min_iter  (int)
%     Minimum number of reweighting iterations.
% param_solver.reweighting_max_iter  (int)
%     Maximum number of reweighting iterations.
% param_solver.reweighting_alpha  (double)
%     Initial value of the reweighting parameter, by default 1.
% param_solver.pdfb_min_iter  (int)
%     Minimum number of iterations.
% param_solver.pdfb_max_iter  (int)
%     Maximum number of iterations.
% param_solver.pdfb_rel_var  (double)
%     Relative variation tolerance.
% param_solver.pdfb_fidelity_tolerance  (double)
%     Tolerance to check data constraints are satisfied.
% param_solver.pdfb_rel_var_low  (double)
%     Minimum relative variation tolerance.
% param_solver.elipse_proj_max_iter  (int)
%     Max. number of iterations, 20.
% param_solver.elipse_proj_min_iter  (int)
%     Min. number of iterations, 1.
% param_solver.elipse_proj_eps  (int)
%     Precision of the projection onto the ellipsoid, 1e-8.
% param_solver.adjust_noise_min_iter (int)
%     Min number of iterations, 100.
% param_solver.adjust_noise_rel_var (double)
%     Tolerance to adjust the noise estimate, 1e-3.
% param_solver.adjust_noise_start_iter (int)
%     Iteration number to force triggering the noise estimation, 500.
% param_solver.adjust_noise_change_percentage (int)
%     The weight of the update w.r.t the l2 norm of the residual data, 0.5.
% param_solver.adjust_noise_start_change_percentage (int)
%     The weight of the update w.r.t the l2 norm of the residual data,
%     if adjustment triggered at ``adjust_noise_start_iter``, 0.1 .
% param_nufft.ox  (int)
%     Oversampling factor (axis x), defaults to 2.
% param_nufft.oy  (int)
%     Oversampling factor (axis y), defaults to 2.
% param_nufft.Kx  (int)
%     Number of neighbour (axis x), defaults to 7.
% param_nufft.Ky  (int)
%     Number of neighbour (axis y), defaults to 7.
% param_nufft.kernel  (string)
%     Name of the nufft interpolation kernel, defaults to "minmax:tuned".
%     Possible options are ``"minmax:tuned"``, ``"minmax:kb"`` or
%     ``"kaiser"``, see ...
% param_precond.N  (int)
%     Number of pixels in the image.
% param_precond.Nox  (int)
%     Number of Fourier points (oversampled plane, x axis).
% param_precond.Noy  (int)
%     Number of Fourier points (oversampled plane, y axis).
% param_precond.gen_uniform_weight_matrix  (bool)
%     Flag to generate uniform weighting matrix (?).
% param_precond.uniform_weight_sub_pixels  (bool)
%     Flag to activate uniform weighting (?).
% dict.wavelet_level  (int)
%     Depth of wavelet decomposition.
% dict.wavelet_basis  (cell{:} of string)
%     Name of the wavelet basis to be used (``"self"`` corresponding to the
%     Dirac basis).
% dict.filter_length (int[:])
%     Length of each wavelet filter considered for the SARA dictionary (0
%     by convention for the Dirac basis).
%
% param_wproj.CEnergyL2 (double)
%     Energy level of the `w`-kernel to be conserved, if `w`-projection is
%     enabled
% param_wproj.GEnergyL2 (double)
%     Energy level  to be conserved of the kernel resulting from the
%     the convolution between the `w`-kernel and the NUFFT kernel, if
%     `w`-projection is enabled.
% param_global_aux.measop_flag_wproj (bool)
%     Flag to activate `w`-correction via `w`-projection, false.
% param_global_aux.measop_flag_visibility_gridding (bool)
%     Flag to activate visibility gridding for data dimensionality
%     reduction, false.
% param_global_aux.adjust_flag_noise (bool)
%     Flag to activate noise estimation on the fly in the case of
%     unreliable noise statistics.
% param_global_aux.parcluster (string)
%     Name of the parallel parcluster profile to launch the parpool. By
%     default  ``"local"`` profile is used. The user should set it to the
%     name of the slurm parcluster profile created on his/her HPC machine,
%     if prefered.
% param_global_aux.algo_flag_computeOperatorNorm (bool)
%     Flag to activate computation of the measurement operator norm,
%     default, true.
% param_global_aux.algo_flag_saveOperatorNorm (bool)
%     Flag to activate saving of the measurement operator norm,
%     default, false.
% param_global_aux.algo_flag_solveMinimization (bool)
%     Flag to trigger the solver, true.
% param_global_aux.preproc_filename_noise_std : anonymous function
%     Function handle, taking the indices of the physical channel
%     and the dataset, and returning a string corresponding to the name of a file
%     containing noise statistics obtained from a pre-processing step.
%     Expected vars: ``sigma``: the standard deviation of the noise,
%     ``RESIDUAL``: (optional) the residual visibilities from a pre-processing
%     step, default  ``[]``.
% param_global_aux.preproc_filename_cal_solutions : anonymous function
%     Function handle, taking the indices of the physical channel and the
%     dataset, and returning a string corresponding to the name of a file
%     containing DDE/DIE calibration solutions.
%     Expected var: ``DDEs``: complex array if DDE calibration and
%     ``DIEs``: complex vector if DIE calibration, default  ``[]``.
% param_global_aux.preproc_filename_model : anonymous function
%     Function handle, taking the indices of the first and last physical channels,
%     associated with an effective (output) channel, and
%     returning the name of a file containing a model image to be used to
%     initialize the reconstruction algorithm, default  ``[]``.
%
% Deprecated fields
%
% param_solver.reweighting_alpha_ff  (double)
%     Update reweighting parameter, if homotopy is active.
%
% param_nnls : struct
%     Structure containing parameters to be passed to the non-negative
%     least-squares (NNLS) algorithm used to estimate the :math:`\ell_2`
%     bounds.
%
% param_nnls.generate_eps_nnls  (bool)
%     Flag to activate NNLS.
% param_nnls.verbose  (int)
%     Print log or not.
% param_nnls.rel_obj  (double)
%     Stopping criterion.
% param_nnls.max_iter  (int)
%     Maximum number of iterations.
% param_nnls.sol_steps  (array of int)
%     Iterations at which the NNLs solution is saved.
% param_nnls.beta  (double)
%     Regularization parameter NNLS.
%
% param_blocking.use_density_partitioning  (bool)
%     Density-based blocking.
% param_blocking.density_partitioning_no  (int)
%     Number of blocks.
% param_blocking.use_uniform_partitioning  (bool)
%     Uniform blocking.
% param_blocking.uniform_partitioning_no  (int)
%     Number of blocks
% param_blocking.use_equal_partitioning  (bool)
%     Equal-size blocking.
% param_blocking.equal_partitioning_no  (int)
%     Number of blocks.
% param_blocking.use_manual_partitioning  (bool)
%     Manual blocking.
% param_blocking.use_manual_frequency_partitioning  (bool)
%     Manual frequency blocking.
%% Parsing json file
fid = fopen(json_filename);
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
config = jsondecode(str);

%% Solver structures
% * Reweighting
% https://fr.mathworks.com/matlabcentral/answers/273955-how-do-i-rename-fields-of-a-structure-array
% if ~isfield(param_global, 'reg_flag_reweighting'); param_global.reg_flag_reweighting = 1; end
% if param_global.reg_flag_reweighting &&  ~isfield(param_global, 'reg_nReweights')
%     param_global.reg_nReweights = 5;
% end
% reweighting_rel_var: double, relative variation (reweighting), 1e-4
% reweighting_min_iter: int,minimum number of reweighting iterations, weights updated reweighting_min_iter times, 5
% reweighting_alpha_ff: double, update reweighting parameter, if homotopy is active. (1 / param_solver.reweighting_alpha)^(1 / (param_solver.reweighting_min_iter - 1));
% reweighting_alpha: double, initial value of the reweighting parameter, 1
param_reweighting = cell2struct(struct2cell(config{1, 1}.reweighting), ...
strcat('reweighting_', fieldnames(config{1, 1}.reweighting)));

% * PDFB
param_pdfb = cell2struct(struct2cell(config{1, 1}.pdfb), ...
strcat('pdfb_', fieldnames(config{1, 1}.pdfb)));

% * Projection onto the ellipsoid
param_proj = cell2struct(struct2cell(config{1, 1}.proj), ...
strcat('elipse_proj_', fieldnames(config{1, 1}.proj)));
% * noise estimation

param_noise = cell2struct(struct2cell(config{1, 1}.noise_estimation), ...
strcat('adjust_noise_', fieldnames(config{1, 1}.noise_estimation)));
% combining the 3 structures (reweighting, pdfb, proj)
param_solver = cell2struct([struct2cell(param_reweighting); struct2cell(param_pdfb); struct2cell(param_proj); struct2cell(param_noise)], ...
    [fieldnames(param_reweighting); fieldnames(param_pdfb); fieldnames(param_proj); fieldnames(param_noise)]);
param_solver.verbose = config{1, 1}.verbose;
param_solver.save_intermediate_results_mat = config{1, 1}.save_intermediate_results_mat;

% * Adaptive epsilon (replaced by noise estimation)
% param_solver.use_adapt_eps = 0;
% % minimum num of iter before stating adjustment
% param_solver.adapt_eps_start = 200;
% % tolerance inside the l2 ball
% param_solver.adapt_eps_tol_in = 0.99;
% % tolerance outside the l2 ball
% param_solver.adapt_eps_tol_out = 1.01;
% % min num of iter between consecutive updates
% param_solver.adapt_eps_steps = 100;
% % bound on the relative change of the solution
% param_solver.adapt_eps_rel_var = 5e-5;
% % the weight of the update w.r.t the l2 norm of the residual data
% param_solver.adapt_eps_change_percentage = (sqrt(5) - 1) / 2;

%% Model structures
% * NUFFT (gridding parameters)
param_nufft = config{2, 1}.nufft;

% * Preconditioning
% set weighting type
param_precond = config{2, 1}.preconditioning; % number of pixels in the image
% ! moved to imaging.m
% % param_precond.N = param_global.im_Nx * param_global.im_Ny; % number of Fourier points (oversampled plane)
% param_precond.Nox = param_nufft.ox * param_global.im_Nx;
% param_precond.Noy = param_nufft.oy * param_global.im_Ny;

% * Wproj
param_wproj = cell2struct(struct2cell(config{2, 1}.wproj), ...
 fieldnames(config{2, 1}.wproj));

% * SARA dictionary
dict = config{2, 1}.sara;

% * Auxilliary global parameters
param_global_aux = cell2struct(struct2cell(config{3, 1}.aux_global), ...
 fieldnames(config{3, 1}.aux_global));
% combining 2 structures (param_global, param_global_aux)
param_global = cell2struct(struct2cell(param_global_aux), ...
     fieldnames(param_global_aux));

end
