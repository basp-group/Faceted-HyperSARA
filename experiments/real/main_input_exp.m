clear; clc; close all;
%% change paths if needed
main_dir =    '/Users/ad33/CodesScience/Faceted-Hyper-SARA/';
project_dir = [main_dir, filesep, 'experiments', filesep, 'real'];
cd(project_dir);
%% src name & datasets
imagecubeName = 'CYGA';
datasetNames = {'CYGA-ConfigA', 'CYGA-ConfigC'}; % allowing for multiple datasets
% data dir.
data_dir = [main_dir, filesep, 'data', filesep, imagecubeName, filesep];
fprintf('\nINFO: data are expected to be saved at %s', data_dir);
% preproc dir.
preproc_calib_dir = [data_dir, 'pre_processing_die/'];
% preproc_calib_dir=[cubedata_dir,'pre_processing_dde/'];

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

idChannels2Image  = [1:16]; % ids of the 'physical' channels to be imaged
nChannelsPerImage   = 16; % number of consecutive channels to be concatenated into each effective channel
nEffChans2Image = floor(numel(idChannels2Image) / nChannelsPerImage); % channels re-structured into effective channels
effChans2Image = cell(nEffChans2Image, 1);
for iEff = 1:nEffChans2Image
    if iEff < nEffChans2Image; effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:iEff * nChannelsPerImage);
    else; effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:end);
    end
    fprintf('\nEffective channel ID %d: physical channels involved: %d - %d', iEff, effChans2Image{iEff}(1), effChans2Image{iEff}(end));
end

%
%% running one subcube at a time
subcubeInd = 0; %  id of subcube if spectral interleaving is active

% measurement op. params
param_global.measop_flag_dataReduction = 1;

% image details, dims &  cellsize
param_global.im_Nx = 2560;
param_global.im_Ny = 1536;
param_global.im_pixelSize = 0.06; % pixelsize in asec

% faceting params: note that if interleaving is active, one subcube is imaged at a time: Qc=1 by default.
param_global.facet_Qx = 1; % dimFacet1
param_global.facet_Qy = 1; % dimFacet2
param_global.facet_overlap_fraction = [0.5, 0.5];

% reg params
param_global.reg_gam = 0.33; % l21 reg param
param_global.reg_gam_bar = 0.33; % nuclear norm reg param

% algo & parallelisation params
param_global.algo_version = 'fhs';

% filenames and input
param_global.exp_type = 'test';
param_global.main_dir = main_dir;
% param_global.preproc_filename_dde = @(firstch,lastch) strcat(calib_dir,'/ddes/chs',num2str(firstch),'-',num2str(lastch),'_dies.mat');
param_global.preproc_filename_die = @(firstch, lastch) strcat(preproc_calib_dir, '/dies/chs', num2str(firstch), '-', num2str(lastch), '_dies.mat');
param_global.preproc_filename_l2bounds = @(firstch, lastch) strcat(preproc_calib_dir, '/l2bounds/chs', num2str(firstch), '-', num2str(lastch), '_l2bounds.mat');
param_global.preproc_filename_model = @(firstch, lastch) strcat(preproc_calib_dir, '/model_images/chs', num2str(firstch), '-', num2str(lastch), '_model_image.fits');

% hardware
param_global.hardware = 'local'; % 'cirrus' or 'local', add your own cluster & update 'util_set_parpool_dev.m' accordingly

%% run main job
main_real_data_exp(imagecubeName, datasetNames, dataFilename, subcubeInd, effChans2Image, param_global);
