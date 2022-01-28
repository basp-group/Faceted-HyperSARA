function [param_global, param_solver, param_nufft, ...
    param_blocking, param_precond, param_nnls, dict] = ...
    read_json_configuration(json_filename, param_global)
% Read algorithmic configuration parameters defined in an input ``.json``
% file.
%
% Parameters
% ----------
% json_filename : string
%     Name of the .json configuration file.
% param_global : struct
%     Global parameter structure containing the following elements.
%
% Returns
% -------
% param_global : struct
%     Structure containing global configuration parameters.
% param_solver : struct
%     Structure containing solver parameters (reweighting scheme, PDFB,
%     projection onto the ellipsoid).
% param_nufft : struct
%     Structure containing the non-uniform FFT parameters.
% param_blocking : struct
%     Structure containing blocking parameters.
% param_precond : struct
%     Structure containing parameters used to define the preconditioning
%     matrix involved in the PDFB algorithm.
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
% param_solver.reweighting_alpha_ff  (double)
%     Update reweighting parameter, if homotopy is active.
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

% Deprecated fields
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

% combining the 3 structures (reweighting, pdfb, proj)
param_solver = cell2struct([struct2cell(param_reweighting); struct2cell(param_pdfb); struct2cell(param_proj)], ...
    [fieldnames(param_reweighting); fieldnames(param_pdfb); fieldnames(param_proj)]);
param_solver.verbose = config{1, 1}.verbose;

if ~isfield(param_global, 'reg_nReweights')
    param_global.reg_nReweights = 5;
end
param_solver.reweighting_max_iter = max(param_global.reg_nReweights, param_solver.reweighting_min_iter + 1);

% * Adaptive epsilon
% TODO : to be removed completely from the solver?
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

%% Model structures
% * NUFFT (gridding parameters)
param_nufft = config{2, 1}.nufft;

% * Blocking
param_blocking = config{2, 1}.blocking;
% partition (symetrically) of the data to nodes (frequency ranges)
param_blocking.fpartition = [icdf('norm', 0.25, 0, pi / 4), 0, ...
    icdf('norm', 0.75, 0, pi / 4), pi];

% * Preconditioning
% set weighting type
param_precond = config{2, 1}.preconditioning;
% number of pixels in the image
param_precond.N = param_global.im_Nx * param_global.im_Ny;
% number of Fourier points (oversampled plane)
param_precond.Nox = param_nufft.ox * param_global.im_Nx;
param_precond.Noy = param_nufft.oy * param_global.im_Ny;

% * NNLS
% if ~isfield(param_global, 'generate_eps_nnls'); param_global.generate_eps_nnls = false; end
% param_nnls = config{2, 1}.nnls;
param_nnls = [];

% * Wproj
% FoV info
% param_wproj.FoVx = sin(pixelSize * Nx * pi / 180 / 3600);
% param_wproj.FoVy = sin(pixelSize * Ny * pi / 180 / 3600);
% param_wproj.uGridSize = 1 / (param_nufft.ox * param_wproj.FoVx);
% param_wproj.vGridSize = 1 / (param_nufft.oy * param_wproj.FoVy);

% * SARA dictionary
% if ~isfield(param_global, 'level'); param_global.wavelet_level = 4;  end
% if ~isfield(param_global, 'basis'); param_global.wavelet_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};  end
dict = config{2, 1}.sara;

end
