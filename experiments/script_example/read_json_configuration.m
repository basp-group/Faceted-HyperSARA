function [param_paths, param_model, speed_of_light, param_solver, param_nufft, param_blocking, param_precond, param_nnls] = read_json_configuration(json_filename)

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
param_paths = config{1, 1};

%% Solver structures
% * Reweighting
param_reweighting = cell2struct(struct2cell(config{3, 1}.reweighting), strcat(fieldnames(config{3, 1}.reweighting), 'reweighting_'));

% * PDFB
% https://fr.mathworks.com/matlabcentral/answers/273955-how-do-i-rename-fields-of-a-structure-array
param_pdfb = cell2struct(struct2cell(config{3, 1}.pdfb), strcat(fieldnames(config{3, 1}.pdfb), 'pdfb_'));

% combining the 2 structures (reweighting and pdfb)
param_solver = [param_pdfb, param_reweighting];

%% Model structures
% * Model
param_model = config{2, 1};

% * NUFFT (gridding parameters)
param_nufft = config{4, 1}.nufft;

% * Blocking
param_blocking = config{4, 1}.blocking;
% partition (symetrically) of the data to nodes (frequency ranges)
param_blocking.fpartition = [icdf('norm', 0.25, 0, pi / 4), 0, icdf('norm', 0.75, 0, pi / 4), pi];

% * Preconditioning
% set weighting type
param_precond = config{4, 1}.preconditioning;
% number of pixels in the image
param_precond.N = Nx * Ny;
% number of Fourier points (oversampled plane)
param_precond.Nox = param_nufft.ox * Nx;
param_precond.Noy = param_nufft.oy * Ny;

% * NNLS
param_nnls = config{4, 1}.nnls;

end
