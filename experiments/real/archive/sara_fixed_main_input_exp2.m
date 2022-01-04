clear; clc; close all; delete(gcp('nocreate'));
%% change path
% cirrus = 1;
% calib_dir = '/lustre/home/shared/sc004/FACETED_HYPERSARA_EXPERIMENTS/calib/SPWIN/'
% cube_filename = '/lustre/home/shared/sc004/FACETED_HYPERSARA_EXPERIMENTS/data/CYG-win-'

cirrus = 1;
main_dir = '//lustre/home/shared/sc004/FACETED_HYPERSARA_EXPERIMENTS/';
calib_dir = [main_dir, 'calib/SPWIN'];
cube_filename = [main_dir, 'data/cyga_data_spwin_'];
project_dir = [main_dir, 'Faceted-Hyper-SARA-EXP2/experiments/real'];

%% general params
% image details, dims &  cellsize
param_global.Nx = 2 * 2560; % floor(2560*sqrt(2)) + mod(floor(2560*sqrt(2)),2);
param_global.Ny = 2 * 1536; % floor(1536*sqrt(2)) + mod(floor(1536*sqrt(2)),2);
param_global.pixelSize = 0.06; % 0.0424;%0.06%0.0424;0.06;%asec
% faceting params
param_global.Qx = 2 * 5;
param_global.Qy = 2 * 3;
param_global.Qc = 1;
param_global.window_type = 'triangular';
param_global.overlap_fraction = [0.5, 0.5];
%
param_global.exp_type = 'sanity';
param_global.main_dir = [main_dir, '/Faceted-Hyper-SARA-EXP2/'];
param_global.die_filename = @(spwin) strcat(calib_dir, '/dies/spwin', num2str(spwin), '_dies.mat');
param_global.l2bounds_filename = @(spwin) strcat(calib_dir, '/l2bounds/spwin', num2str(spwin), '_l2bounds.mat');
param_global.model_filename = @(spwin) strcat(calib_dir, '/model_image/spwin', num2str(spwin), '_model_image.fits');

%% flags
input_flags.computeOperatorNorm = 1;
input_flags.solveMinimization = 1;
input_flags.homotopy = 0;
input_flags.dr = 1;
input_flags.wprojection = 0;

% cd ([project_dir,filesep,'experiments',filesep,'real']);
imagename = 'CygA';

spwins2image = 31; [1:17 19 21:32]; % ids of the channels to be imaged from spwin=spwin_ind;
% if 'sara', make sure to select the id of the channel to be imaged.

ncoredata =  numel(spwins2image); % 30 for cyga

algoversion = 'sara';
subcube_ind = 1;
param_reg.gam = 1;
param_reg.gam_bar = 1;
param_reg.rw = 1;
param_reg.nReweights = 5;

main_real_data_dev_exp22(imagename, cube_filename, subcube_ind, spwins2image, ...
    algoversion, ncoredata, param_global, param_reg, input_flags, cirrus);
