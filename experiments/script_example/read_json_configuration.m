function [speed_of_light, param_general, param_model, param_solver, ...
    param_nufft, param_blocking, param_precond, param_synth, param_nnls] = ...
    read_json_configuration(json_filename, ind)

% TODO: modify the name of the model parameters (use a structure to define 
% these)

%% Parsing json file
fid = fopen(json_filename);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
config = jsondecode(str);

speed_of_light = 299792458;

%% Paths (to dataset and uv-coverage)
param_general = config{1, 1}.general;

%% (Optional) Parameters for synthetic data generation (only useful for 
% synthetic data scirpts)
param_synth = config{1, 1}.synth;

%% Solver structures
% * Reweighting
% https://fr.mathworks.com/matlabcentral/answers/273955-how-do-i-rename-fields-of-a-structure-array
param_reweighting = cell2struct(struct2cell(config{2, 1}.reweighting), strcat('reweighting_', fieldnames(config{2, 1}.reweighting)));

% * PDFB
param_pdfb = cell2struct(struct2cell(config{2, 1}.pdfb), strcat('pdfb_', fieldnames(config{2, 1}.pdfb)));

% * Projection onto the ellipsoid
param_proj = cell2struct(struct2cell(config{2, 1}.proj), strcat('elipse_proj_', fieldnames(config{2, 1}.proj)));

% combining the 2 structures (reweighting and pdfb)
param_solver = cell2struct([struct2cell(param_reweighting); struct2cell(param_pdfb);struct2cell(param_proj)], ...
    [fieldnames(param_reweighting);fieldnames(param_pdfb);fieldnames(param_proj)]);
param_solver.cube_id = ind;

% * Adaptive epsilon
% TODO : to be removed form the solver
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
% * Model
param_model = config{1, 1}.mandatory;

% * NUFFT (gridding parameters)
param_nufft = config{3, 1}.nufft;

% * Blocking
param_blocking = config{3, 1}.blocking;
% partition (symetrically) of the data to nodes (frequency ranges)
param_blocking.fpartition = [icdf('norm', 0.25, 0, pi / 4), 0, icdf('norm', 0.75, 0, pi / 4), pi];

% * Preconditioning
% set weighting type
param_precond = config{3, 1}.preconditioning;
% number of pixels in the image
param_precond.N = param_model.Nx * param_model.Ny;
% number of Fourier points (oversampled plane)
param_precond.Nox = param_nufft.ox * param_model.Nx;
param_precond.Noy = param_nufft.oy * param_model.Ny;

% * NNLS
param_nnls = config{3, 1}.nnls;

end
