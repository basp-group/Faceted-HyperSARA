clear; clc; close all;
%% change paths if needed
% TODO (keep empty parameter for non-used variables)
main_dir = '../..'; % '/Users/ad33/CodesScience/Faceted-Hyper-SARA/';
project_dir = [main_dir, filesep, 'experiments', filesep, 'spatial'];
cd(project_dir);
%% src name & datasets
imagecubeName = 'cygASband_Cube_256_512';
param_global.exp_type = 'test';
datasetNames = {'test'}; % allowing for multiple datasets
% data dir.
data_dir = [main_dir, filesep, 'data', filesep, imagecubeName, filesep];
fprintf('\nINFO: data are expected to be saved at %s\n', data_dir);
% preproc dir.
preproc_calib_dir = [data_dir, 'pre_processing_die/'];
% preproc_calib_dir=[cubedata_dir,'pre_processing_dde/'];
% name of json parameter file
json_filename = "default_parameters.json";

% data files
% example of data file: 'data_ch_1.mat' with vars : 'maxProjBaseline', 'u','v','w', 'nW', 'y', 'frequency'.
% Note that data 'y' are not whitened, uvw coordinates are in units of the
% wavelength (i.e. normalised with the freq) and 'maxProjBaseline' is in
% units of the wavelength
% ! Note the path for .mat files (dataset nametag used)
dataFilename = @(idSet, ch) strcat(data_dir, filesep, datasetNames{idSet}, filesep, 'data_ch_', num2str(ch), '.mat');

%% channels organisation
%%%% option 1: provide a cell array containing the ids of the  channels to be concatenated for each effective channel.
% example a: two effective channels, containing two 'physical' channels each
% > effChans2Image={[1,2],[3,4]};

% example b: one channel effective channel with one physical channel
% > effChans2Image={[1]}

%%%% option 2: provide all ids of channels 'nChannelsPerImage' & num of channel per effective channel 'nChannelsPerImage' channel
% example c: EXP 1: subcube 1 (first channel from each spectral window (2 of  which are flagged).
% >idChannels2Image  = [1:16:272 289:16:320 337:16:512];
% >nChannelsPerImage   = 1;

% example d: Exp 2: reduced imagecube containing 30 effective channels each concatenating 16 physical channels
% >idChannels2Image  = [1:272 289:320 337:512];
% >nChannelsPerImage   = 16;

% TODO
idChannels2Image  = [1:2]; % ids of the 'physical' channels to be imaged
nChannelsPerImage   = 1; % number of consecutive channels to be concatenated into each effective channel
nEffChans2Image = floor(numel(idChannels2Image) / nChannelsPerImage); % channels re-structured into effective channels
effChans2Image = cell(nEffChans2Image, 1);
for iEff = 1:nEffChans2Image
    if iEff < nEffChans2Image; effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:iEff * nChannelsPerImage);
    else; effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:end);
    end
    fprintf('\nEffective channel ID %d: physical channels involved: %d - %d\n', iEff, effChans2Image{iEff}(1), effChans2Image{iEff}(end));
end

%
%% running one subcube at a time
% TODO
subcubeInd = 0; %  id of subcube if spectral interleaving is active, ([P.A.] 0 otherwise?)

% measurement op. params
param_global.measop_flag_dataReduction = 0;

% image details, dims &  cellsize
% TODO
param_global.im_Nx = 512; % 2560;
param_global.im_Ny = 256; % 1536;
param_global.im_pixelSize = []; % 0.06; % pixelsize in asec ([P.-A.] = [] if using synth data)
param_global.algo_flag_computeOperatorNorm = false;

% faceting params: note that if interleaving is active, one subcube is imaged at a time: Qc=1 by default.
param_global.facet_Qx = 1; % dimFacet1
param_global.facet_Qy = 1; % dimFacet2
param_global.facet_overlap_fraction = [0.5, 0.5];

% reg params
param_global.reg_gam = 0.33; % l21 reg param
param_global.reg_gam_bar = 0.33; % nuclear norm reg param
param_global.reg_flag_reweighting = true;
param_global.reg_flag_homotopy = false;
param_global.reg_nReweights = 5;
param_global.generate_eps_nnls = false;

% algo & parallelisation params
param_global.algo_version = 'fhs'; % 'fhs', 'hs' or 'sara'

% filenames and input
param_global.main_dir = main_dir;
% param_global.preproc_filename_dde = @(firstch,lastch) strcat(calib_dir,'/ddes/chs',num2str(firstch),'-',num2str(lastch),'_dies.mat');
% param_global.preproc_filename_die = @(firstch, lastch) strcat(preproc_calib_dir, '/dies/chs', num2str(firstch), '-', num2str(lastch), '_dies.mat');
param_global.preproc_filename_die = [];
% param_global.preproc_filename_l2bounds = @(firstch, lastch) strcat(preproc_calib_dir, '/l2bounds/chs', num2str(firstch), '-', num2str(lastch), '_l2bounds.mat');
param_global.preproc_filename_l2bounds = [];
% param_global.preproc_filename_model = @(firstch, lastch) strcat(preproc_calib_dir, '/model_images/chs', num2str(firstch), '-', num2str(lastch), '_model_image.fits');
param_global.preproc_filename_model = [];

% hardware
param_global.hardware = 'local'; % 'cirrus' or 'local', add your own cluster & update 'util_set_parpool_dev.m' accordingly

%% read and set configuration from .json file
[param_global, param_solver, ...
    param_nufft, param_blocking, param_precond, param_nnls, dict] = ...
    read_json_configuration(json_filename, param_global);

%% run main job
main_real_data_exp(imagecubeName, datasetNames, dataFilename, ...
    subcubeInd, effChans2Image, param_solver, param_global, ...
    param_nufft, param_blocking, param_precond, param_nnls, dict);
