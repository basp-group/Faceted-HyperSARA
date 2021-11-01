function [speed_of_light, param_paths, param_model, param_solver, ...
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
param_paths = config{1, 1}.general;

%% (Optional) Parameters for synthetic data generation (only useful for 
% synthetic data scirpts)
param_synth = config{1, 1}.synth;

%% Solver structures
% * Reweighting
% https://fr.mathworks.com/matlabcentral/answers/273955-how-do-i-rename-fields-of-a-structure-array
param_reweighting = cell2struct(struct2cell(config{3, 1}.reweighting), strcat(fieldnames(config{2, 1}.reweighting), 'reweighting_'));

% * PDFB
param_pdfb = cell2struct(struct2cell(config{2, 1}.pdfb), strcat(fieldnames(config{2, 1}.pdfb), 'pdfb_'));

% * Projection onto the ellipsoid
param_proj = cell2struct(struct2cell(config{2, 1}.proj), strcat(fieldnames(config{2, 1}.proj), 'elipse_proj_'));

% combining the 2 structures (reweighting and pdfb)
param_solver = [param_reweighting, param_pdfb, param_proj];
param_solver.cube_id = ind;

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
param_precond.N = Nx * Ny;
% number of Fourier points (oversampled plane)
param_precond.Nox = param_nufft.ox * Nx;
param_precond.Noy = param_nufft.oy * Ny;

% * NNLS
param_nnls = config{3, 1}.nnls;

end
