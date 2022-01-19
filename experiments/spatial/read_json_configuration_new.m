function [speed_of_light, param_global, param_solver, ...
    param_nufft, param_blocking, param_precond, param_nnls, dict] = ...
    read_json_configuration_new(json_filename, param_global)

% TODO: modify the name of the model parameters (use a structure to define
% these)
% TODO: include full list of fields set in this function

%% Parsing json file
fid = fopen(json_filename);
raw = fread(fid, inf);
str = char(raw');
fclose(fid);
config = jsondecode(str);

speed_of_light = 299792458;

%% Solver structures
% * Reweighting
% https://fr.mathworks.com/matlabcentral/answers/273955-how-do-i-rename-fields-of-a-structure-array
% ! revise place of parameter
% if ~isfield(param_global, 'reg_flag_reweighting'); param_global.reg_flag_reweighting = 1; end
% if param_global.reg_flag_reweighting &&  ~isfield(param_global, 'reg_nReweights')
%     param_global.reg_nReweights = 5;
% end
% reweighting_rel_var: double, relative variation (reweighting), 1e-4
% reweighting_min_iter: int,minimum number of reweighting iterations, weights updated reweighting_min_iter times, 5
% reweighting_alpha_ff: double, update reweighting parameter, if homotopy is active. (1 / param_solver.reweighting_alpha)^(1 / (param_solver.reweighting_min_iter - 1));
% reweighting_alpha: double, initial value of the reweighting parameter, 1
param_reweighting = cell2struct(struct2cell(config{1, 1}.reweighting), strcat('reweighting_', fieldnames(config{1, 1}.reweighting)));

% * PDFB
% pdfb_min_iter: int, minimum number of iterations
% pdfb_max_iter: int, maximum number of iterations
% pdfb_rel_var: double, relative variation tolerance
% pdfb_fidelity_tolerance: double, tolerance to check data constraints are satisfied
% pdfb_rel_var_low: double, minimum relative variation tolerance 
param_pdfb = cell2struct(struct2cell(config{1, 1}.pdfb), strcat('pdfb_', fieldnames(config{1, 1}.pdfb)));

% * Projection onto the ellipsoid
% elipse_proj_max_iter: int, max. number of iterations, 20
% elipse_proj_min_iter: int, min. number of iterations, 1
% elipse_proj_eps: int, precision of the projection onto the ellipsoid, 1e-8
param_proj = cell2struct(struct2cell(config{1, 1}.proj), strcat('elipse_proj_', fieldnames(config{1, 1}.proj)));

% combining the 3 structures (reweighting, pdfb, proj)
param_solver = cell2struct([struct2cell(param_reweighting); struct2cell(param_pdfb); struct2cell(param_proj)], ...
    [fieldnames(param_reweighting); fieldnames(param_pdfb); fieldnames(param_proj)]);
param_solver.verbose = config{1, 1}.verbose;

if param_global.reg_flag_reweighting &&  ~isfield(param_global, 'reg_nReweights')
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
% ox: int, oversampling factor (x), defaults to 2
% oy: int, oversampling factor (y), defaults to 2
% Kx: int, number of neighbour (x), defaults to 7
% Ky: int, number of neighbour (y), defaults to 7
% kernel: string, nufft interpolation kernel, defaults to 'minmax:tuned'
param_nufft = config{2, 1}.nufft;

% * Blocking
% % density-based
% use_density_partitioning
% density_partitioning_no
% % uniform
% use_uniform_partitioning
% uniform_partitioning_no
% % equal-size
% use_equal_partitioning
% equal_partitioning_no
% % manual
% use_manual_partitioning
% use_manual_frequency_partitioning
% % partition (symetrically) of the data to nodes (frequency ranges)
% fpartition = [icdf('norm', 0.25, 0, pi / 4), 0, icdf('norm', 0.75, 0, pi / 4), pi];
%
param_blocking = config{2, 1}.blocking;
% partition (symetrically) of the data to nodes (frequency ranges)
param_blocking.fpartition = [icdf('norm', 0.25, 0, pi / 4), 0, icdf('norm', 0.75, 0, pi / 4), pi];

% * Preconditioning
% N: number of pixels in the image
% Nox: number of Fourier points (oversampled plane)
% Noy:number of Fourier points (oversampled plane)
% set weighting type
% gen_uniform_weight_matrix
% uniform_weight_sub_pixels
%
% set weighting type
param_precond = config{2, 1}.preconditioning;
% number of pixels in the image
param_precond.N = param_global.im_Nx * param_global.im_Ny;
% number of Fourier points (oversampled plane)
param_precond.Nox = param_nufft.ox * param_global.im_Nx;
param_precond.Noy = param_nufft.oy * param_global.im_Ny;

% * NNLS
% generate_eps_nnls: flag to activate NNLS
% verbose: print log or not
% rel_obj: stopping criterion
% max_iter: max number of iterations
% sol_steps: (vector) saves images at the given iterations
% beta: double, regularization parameter 
% ! revise place of parameter
% if ~isfield(param_global, 'generate_eps_nnls'); param_global.generate_eps_nnls = false; end
param_nnls = config{2, 1}.nnls;

% * Wproj
% FoV info
% param_wproj.FoVx = sin(pixelSize * Nx * pi / 180 / 3600);
% param_wproj.FoVy = sin(pixelSize * Ny * pi / 180 / 3600);
% param_wproj.uGridSize = 1 / (param_nufft.ox * param_wproj.FoVx);
% param_wproj.vGridSize = 1 / (param_nufft.oy * param_wproj.FoVy);

% * SARA dictionary
% ! revise place of parameter
% if ~isfield(param_global, 'level'); param_global.wavelet_level = 4;  end
% if ~isfield(param_global, 'basis'); param_global.wavelet_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};  end
dict = config{2, 1}.sara;

end
