clear; clc; close all; delete(gcp('nocreate'));

%% change path
cirrus = 0;
main_dir = '/Users/ad33/CodesScience/FHS_EXP12/Faceted-Hyper-SARA/';
calib_dir = [main_dir, 'calib'];
cube_filename = [main_dir, 'data/CYG/cyga_data_subcube_'];
project_dir = [main_dir, 'experiments/real'];

%% general params
% image details, dims &  cellsize
param_global.Nx = 2560;
param_global.Ny = 1536;
param_global.pixelSize = 0.06; % asec
% faceting params
param_global.Qx = 1; 2 * 5;
param_global.Qy = 1; 2 * 3;
param_global.Qc = 1;
param_global.window_type = "triangular";
param_global.overlap_fraction = [0.25, 0.25];
%
param_global.exp_type = 'realexp';
param_global.main_dir = main_dir;
param_global.G_filename = @(subcube, ch) strcat(calib_dir, '/Gw/subcube', num2str(subcube), '_ch', num2str(ch), '_Gw.mat');
param_global.l2bounds_filename = @(subcube, ch) strcat(calib_dir, '/l2bounds/subcube', num2str(subcube), '_ch', num2str(ch), '_l2bounds.mat');
param_global.model_filename = @(subcube, ch) strcat(calib_dir, '/model_image/subcube', num2str(subcube), '_ch', num2str(ch), '_model_image.fits');

%% flags
input_flags.computeOperatorNorm = 1;
input_flags.solveMinimization = 1;
input_flags.dr = 0;
input_flags.wprojection = 0;

%% params

channels2image = 1; [1:17 19 21:32]; % ids of the channels to be imaged from subcube=subcube_ind;
% if 'sara', make sure to select the id of the channel to be imaged.

ncoredata = numel(channels2image); % 30 for cyga
ncoredata = param_global.Qx * param_global.Qy + numel(channels2image) + 1; % 30 for cyga

fprintf('\nINFO:  %d cores deployed !\n', ncoredata);

%% local run
if ~cirrus
    imagename = 'CYG';
    algoversion = 'sara';
    subcube_ind = 1;
    param_reg.gam = 1;
    param_reg.gam_bar = 1;
    param_reg.rw = 1;
    param_reg.nReweights = 5;
    % main script"  >>$RUNFILENAME
    main_real_data_dev_exp1(imagename, cube_filename, subcube_ind, channels2image, algoversion, ncoredata, param_global, param_reg, input_flags, cirrus);
end
