% Main imaging utility script, to be configured by the user to run an
% experiment.
%
% This file calls the generic imaging pipeline
% :mat:func:`imaging.imaging` to reconstruct an image from
% the input dataset.
%
% Parameters
% ----------
% srcName : string
%     Name of the target source to be reconstructed as specified during
%     data extraction. It will be used for output files.
% datasetNames : cell of string
%     Name of the different datasets used (e.g.,
%     different acquisition configurations, different observation times),
%     to be set to ``{''}`` if one data set is imaged or no name tag of 
%     the MS is given during data extraction.
% json_filename : string
%     Name of the input ``.json`` configuration file specifying the
%     value of the algorithm parameters (PDFB, reweighting, ellipsoid,
%     projection, preconditioning, ...).
% dataFilename : anonymous function
%     Function handle taking a channel index as an input, and returning the
%     name of the associated data ``.mat`` file (one file per wavelength), nomenclature adopted during data extraction .
% idChannels2Image : array[int]
%     Indices of the `physical` channels to be imaged.
% nChannelsPerImage : int
%     Number of consecutive channels to be concatenated into each effective
%     channel (1 if full spectral resolution is considered, i.e. number of output channels is equal to the number of input channels).
% subcubeInd : int
%     Index of the subcube to be reconstructed (if spectral faceting is
%     used). To be set to ``0`` otherwise.
% main_dir : string
%     Directory of the Faceted-HyperSARA repository.
% data_dir : string
%     Directory where the input ``.mat`` data files are saved.
% project_dir : string
%     Directory of the experiment project.
% preproc_calib_dir : string
%     Name of the folder containing results from a calibration
%     pre-processing step. To be set to ``[]`` if not used or not
%     available.
% param_global.preproc_filename_die : anonymous function
%     Function handle, taking the index of the first and last channel, and
%     returning a string corresponding to the name of a file containing
%     DIE calibration constants. If not used, must be set to ``[]``.
% param_global.preproc_filename_l2bounds : anonymous function
%     Function handle, taking the index of the first and last channel, and
%     returning a string corresponding to the name of a file
%     containing pre-computed :math:`\ell_2` bounds. If not used, must be set
%     to ``[]``.
% param_global.preproc_filename_model : anonymous function
%     Function handle, taking the index of the first and last channel, and
%     returning the name of a file containing a model image to be used to
%     initialize the reconstruction algorithm. If not used, should be set to
%     ``[]``.
% param_global.measop_flag_visibility_gridding : bool
%     Flag to activate visibility gridding for data dimensionality reduction
%     :cite:p:`Kartik2017`.
% param_global.algo_flag_computeOperatorNorm : bool
%     Flag to trigger the computation of the norm of the measurement
%     operator. If set to ``false``, MATLAB will look for a file where this
%     quantity has been saved (save and computation step are triggered in
%     :mat:func:`imaging.imaging`).
% param_global.im_Nx  : int
%     Number of pixels along axis x in the reconstructed image.
% param_global.im_Ny : int
%     Number of pixels along axis y in the reconstructed image.
% param_global.im_pixelSize  : double
%     Pixel-size in arcsec. Set to ``[]`` to use
%     the default computation, that is 2x the resolution of the
%     observation.
% param_global.reg_flag_reweighting : bool
%     Flag to activate the re-weighting procedure.
% param_global.reg_nReweights : int
%     Maximum number of reweighting steps.
% param_global.parcluster : string
%     String the name of the parallel parcluster profile. By default ``"local"`` profile
%     is used. The user should set it to the name of the slurm parcluster
%     profile created on his/her HPC machine, if prefered. 
% param_global.algo_solver : string
%     Name of the imaging approach used. Possible values are ``"fhs"``
%     (Faceted-HyperSARA), ``"hs"`` (HyperSARA) and ``"sara"`` (SARA).
% param_global.facet_Qx : int
%     Number of spatial facets along spatial axis x. Will be reset
%     automatically to 1 if ``param_global.algo_solver = "sara"`` or
%     ``"hs"``.
% param_global.facet_Qy : int
%     Number of spatial facets along spatial axis y. Will be reset
%     automatically to 1 if ``param_global.algo_solver = "sara"`` or
%     ``"hs"``.
% param_global.facet_overlap_fraction : array[double]
%     Array containing the overlap fraction between consecutive facets along
%     each axis (y and x) for the faceted low-rankness prior. Will be reset
%     automatically to ``[0, 0]`` if ``param_global.algo_solver = "sara"``
%     or ``"hs"``. Besides, each entry of
%     ``param_global.facet_overlap_fraction`` is reset to 0 the
%     number of facet along the corresponding dimension is equal to 1
%     (i.e., ``param_global.facet_Qy = 1`` and/or
%     ``param_global.facet_Qx = 1``).
% param_global.reg_gam : double
%     Average joint sparsity prior regularization parameter.
% param_global.reg_gam_bar : double
%     Low-rankness prior regularization parameter (for HyperSARA and
%     Faceted HyperSARA).
%
% Example
% -------
%
% .. code-block:: matlab
%
%    %% Set the input and output channels IDs
%    %% option 1: provide a cell array containing the ids of the channels
%    % to be concatenated for each effective channel.
%    % example a: two effective channels, containing two "physical"
%    % channels each
%    effChans2Image = {[1,2], [3,4]};
%
%    % example b: one channel effective channel with one physical channel
%    effChans2Image = {[1]};
%
%
%    %% option 2: provide all ids of channels nChannelsPerImage & number
%    % of channel per effective channel nChannelsPerImage channel
%    % example c: EXP 1: subcube 1 (first channel from each spectral
%    % window (2 of  which are flagged).
%    idChannels2Image = [1:16:272, 289:16:320, 337:16:512];
%    nChannelsPerImage = 1;
%
%    % example d: reduced imagecube containing 30 effective channels
%    % each concatenating 16 physical channels
%    idChannels2Image = [1:272, 289:320, 337:512];
%    nChannelsPerImage = 16;
%
% Warning
% -------
% - Data files (example: ``'data_ch_1.mat'``) are expected to contain the following variables : `maxProjBaseline`
%   (double, maximum projected baseline), `u` , `v`, `w` (arrays of
%   :math:`uvw`-coordinates), `nW` (noise whitening vector), `y` (complex vector of
%   visibilities, Stokes I), `frequency` (associated frequency value).
%
% - Note that the data `y`  are not whitened. 
%   The variables `maxProjBaseline` and the `uvw` coordinates should be 
%   in units of the wavelength (i.e. normalised by the associated wavelength).
%

%% Documentation
%
% For more details about the different parameters to be configured, see the
% online documentation: https://basp-group.github.io/Faceted-Hyper-SARA/default.html .

clear; clc; close all;

%% Parameters to be specified
% TODO (keep empty parameter for non-used variables)

% change paths if needed
% TODO: to be adjusted by the user if required (keep current setting by default)
main_dir = '.'; % Faceted-Hyper-SARA dir.
project_dir = [main_dir, filesep, 'imaging'];
cd(project_dir);

% src name & datasets
% TODO: to be adjusted by the user
srcName = 'CYGA'; % as specified during data extraction
datasetNames = {'CYGA-ConfigA', 'CYGA-ConfigC'}; %  accomodating multiple datasets, 
%set as specified in the data extraction step, {''} otherwise.

% data directory (as expected from the data extraction step)
data_dir = [main_dir, filesep, 'data', filesep, srcName, filesep];
fprintf('\nINFO: data are expected to be saved at %s\n', data_dir);

% data files
% example of data file: `./data/CYGA/CYGA-ConfigA/data_ch_1.mat` 
% ! Note in the path for .mat files (dataset name tag used)
% ! To be adjusted only if data paths during data extraction have been altered
dataFilename = @(idSet, ch) strcat(data_dir, filesep, datasetNames{idSet}, filesep, 'data_ch_', num2str(ch), '.mat');


% calibration preprocessing step directory
% TODO: to be adjusted by the user (set path to [] whenever it is not used)
preproc_calib_dir = [data_dir,filesep, 'pre_processing_calibration/']; 
% check if directory of calibration results exists
if exist(preproc_calib_dir, 'dir')
    preproc_results_exist = 1;
else,  preproc_results_exist = 0;
end
% name of json parameter file
json_filename = "default_parameters.json";

%% channels organisation
%%%% option 1: provide a cell array containing the ids of the  channels to be concatenated for each effective channel.
% example a: two effective channels, containing two 'physical' channels each
% > effChans2Image={[1,2],[3,4]};

% example b: one channel effective channel with one physical channel
% > effChans2Image={[1]}

%%%% option 2: provide all ids of channels 'nChannelsPerImage' & num of channel per effective channel 'nChannelsPerImage' channel
% 'effChans2Image' will then be determined automatically.
% example c: EXP 1: subcube 1 (first channel from each spectral window (2 of  which are flagged).
% > idChannels2Image = [1:16:272 289:16:320 337:16:512];
% > nChannelsPerImage = 1;

% example d: Exp 2: reduced imagecube containing 30 effective channels each concatenating 16 physical channels
% > idChannels2Image = [1:272 289:320 337:512];
% > nChannelsPerImage = 16;

% TODO: to be adjusted by the user
idChannels2Image = [1:256]; % ids of the 'physical' channels to be imaged, e.g., [1:256] for channel range from 1 to 256
nChannelsPerImage = 16; % number of consecutive channels to be combined into each effective channel

nEffChans2Image = floor(numel(idChannels2Image) / nChannelsPerImage); % ouput effective channels: number of images in the estimate model cube
% !must be equal to 1 for 'sara' and greater than 1 for 'fhs'

%re-arranging channels indices based on idChannels2Image, nChannelsPerImage and nEffChans2Image
effChans2Image = cell(nEffChans2Image, 1);
for iEff = 1:nEffChans2Image
    if iEff < nEffChans2Image; effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:iEff * nChannelsPerImage);
    else; effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:end);
    end
    fprintf('\nINFO:Effective channel ID %d: physical channels involved: %d - %d', iEff, effChans2Image{iEff}(1), effChans2Image{iEff}(end));
end

% running one subcube at a time
% TODO: to be adjusted by the user
subcubeInd = 0; % id of subcube if spectral interleaving is active, 0 if inactive

% measurement operator features
% TODO: to be adjusted by the user
% activating visibility gridding to reduce data dimensionality
param_global.measop_flag_visibility_gridding = 0; % 1 if active, 0 otherwise. 

% image details, dims &  cellsize
% TODO: to be adjusted by the user
param_global.im_Nx = 2560;
param_global.im_Ny = 1536;
param_global.im_pixelSize =  0.06; % pixelsize in asec, set to [] to use the default value set from the uv-coverage

% faceting params: note that if interleaving is active, one subcube is
% imaged at a time: Qc=1 by default.
% TODO: to be adjusted by the user
param_global.facet_Qx = 5; % number of spatial facets along axis x (only used in Faceted HyperSARA)
param_global.facet_Qy = 3; % number of spatial facets along axis y (only used in Faceted HyperSARA)
param_global.facet_overlap_fraction = [0.1, 0.1];

% reg params
% TODO: to be adjusted by the user
param_global.reg_gam = 0.33; % additional scaling factor for sparsity 
% regularization (using heuristics described in Thouvenin2021). Set parameter to 1 by default.
param_global.reg_gam_bar = 0.33; % additional scaling factor for low-rankness regularization (using heuristics described in Thouvenin2021). Only active 
% for HyperSARA and Faceted HyperSARA. Set parameter to 1 by default.
param_global.reg_flag_reweighting = true;
param_global.reg_nReweights = 5;

% algo & parallelisation params
% TODO: to be adjusted by the user
param_global.algo_solver = 'fhs'; % 'fhs' (Faceted HyperSARA), 'hs' (HyperSARA) or 'sara' (SARA approach)

% filenames and input
param_global.main_dir = main_dir;

if preproc_results_exist
    % TODO: to be adjusted by the user (set path to [] whenever it is not used)
    param_global.preproc_filename_dde = @(firstch, lastch) strcat(preproc_calib_dir, filesep, 'ddes', filesep, 'chs', num2str(firstch), '-', num2str(lastch), '_ddes.mat');
    param_global.preproc_filename_die = @(firstch, lastch) strcat(preproc_calib_dir, filesep,'dies',filesep,'chs', num2str(firstch), '-', num2str(lastch), '_dies.mat');
    param_global.preproc_filename_l2bounds = @(firstch, lastch) strcat(preproc_calib_dir, filesep, 'l2bounds', filesep, 'chs', num2str(firstch), '-', num2str(lastch), '_l2bounds.mat');
    param_global.preproc_filename_model = @(firstch, lastch) strcat(preproc_calib_dir, filesep, 'model_images', filesep, 'chs', num2str(firstch), '-', num2str(lastch), '_model_image.fits');
end


% Parallel parcluster profile 
% TODO: to be adjusted by the user
param_global.parcluster = 'local'; % name of the parallel parcluster profile to use.
% to be adjusted if running on a HPC and a slurm parcluster profile will be used. 
%% read and set configuration from .json file
[param_global, param_solver, ...
    param_nufft, param_precond, dict] = ...
    read_json_configuration(json_filename, param_global);

%% run main job
imaging(srcName, datasetNames, dataFilename, ...
    subcubeInd, effChans2Image, param_solver, param_global, ...
    param_nufft, param_precond, dict);
