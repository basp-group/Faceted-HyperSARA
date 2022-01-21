% Main imaging utility script, to be configured by the user to run an
% experiment.
%
% This file calls the generic imaging pipeline
% :func:`experiments.main_real_data_exp` to reconstruct an image from the
% input dataset.
%
% Parameters
% ----------
% imagecubeName : string
%     Name of the image cube to be reconstructed.
% datasetNames : cell of string
%     Name of the different datasets used (e.g., if datasets associated with
%     different acquisition configurations are used.
% json_filename : string
%     Name of the input ``.json`` configuration file specifying the
%     value of the algorithm parameters (PDFB, reweighting, ellipsoid,
%     projection, preconditioning, ...).
% dataFilename : anonymous function
%     Function handle taking a channel index as an input, and returning the
%     name of the associated data ``.mat`` file (one file per wavelength).
% idChannels2Image : array[int]
%     Indices of the `physical` channels to be imaged.
% nChannelsPerImage : int
%     Number of consecutive channels to be concatenated into each effective
%     channel (1 if no reduction is used).
% subcubeInd : int
%     Index of the subcube to be reconstructed (if spectral faceting is
%     used). To be set to ``0`` otherwise.
% main_dir : string
%     Directory of the Faceted-HyperSARA repository.
% data_dir : string
%     Directory where the input ``.mat`` data files are saved.
% project_dir : string
%     Directory of the experiment project.
% preproc_calib_dir : anonymous function
%     Name of the folder containing pre-estimated calibration kernels. If
%     not used or not available, to be set to ``[]``
% param_global.preproc_filename_die : anonymous function
%     Function handle, taking the index of the first and last channel, and
%     returning a string corresponding to the name of a file containing
%     DIE calibration constants. If not used, to be set to ``[]``.
% param_global.preproc_filename_l2bounds : anonymous function
%     Function handle, taking the index of the first and last channel, and
%     returning a string corresponding to the name of a file
%     containingpre-computed :math:`ell_2` bounds. If not used, to be set
%     to ``[]``.
% param_global.preproc_filename_model : anonymous function
%     Function handle, taking the index of the first and last channel, and
%     returning the name of a file containing a model image to be used to
%     initialize the reconstruction algorithm. If not used, to be set to
%     ``[]``.
% param_global.measop_flag_dataReduction : bool
%     Flag to activate data reduction.
% param_global.algo_flag_computeOperatorNorm : bool
%     Flag to trigger the computation of the norm of the measurement
%     operator. If set to ``false``, MATLAB will look for a file where this
%     quantity has been saved (save and computation step are triggered in
%     :func:`main_real_data_exp`).
% param_global.im_Nx  : int
%     Number of pixels along axis x in the reconstructed image.
% param_global.im_Ny : int
%     Number of pixels along axis y in the reconstructed image.
% param_global.im_pixelSize  : double
%     Pixel-size in arcsec in the frequency domain. Set to ``[]`` to use
%     the default computation (implemented in :func:`main_real_data_exp`).
% param_global.generate_eps_nnls : bool
%     Flag to activate the computation of the :math:`\ell_2` bounds.
% param_global.reg_flag_reweighting : bool
%     Flag to activate the use of multiple reweighting steps.
% param_global.reg_nReweights : int
%     Maximum number of reweighting steps.
% param_global.reg_flag_homotopy : bool
%     Flag to use the reweighting homotopy strategy (deactivated for now).
% param_global.hardware : string
%     String to specify the hardware configuration for the parallelization.
%     Possible values are ``'cirrus'`` or ``'local'``.
% param_global.algo_version : string
%     Name of the imaging approach used. Possible values are ``'fhs'``
%     (Faceted-HyperSARA), ``'hs'`` (HyperSARA) and ``'sara'`` (SARA).
% param_global.facet_Qx : int
%     Number of spatial facets along spatial axis x. Will be reset
%     automaticallt to 1 if ``param_global.algo_version = 'sara'`` or 'hs'.
% param_global.facet_Qy : int
%     Number of spatial facets along spatial axis y. Will be reset
%     automaticallt to 1 if ``param_global.algo_version = 'sara'`` or 'hs'.
% param_global.facet_overlap_fraction : array[double]
%     Array containing the overlap fraction between consecutive facets along
%     each axis (y and x). Will be reset
%     automatically to ``[0, 0]`` if ``param_global.algo_version = 'sara'``
%     or 'hs'. Besides, each entry of
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
%    %% option 1: provide a cell array containing the ids of the channels
%    % to be concatenated for each effective channel.
%    % example a: two effective channels, containing two 'physical'
%    % channels each
%    effChans2Image = {[1,2],[3,4]};
%
%    % example b: one channel effective channel with one physical channel
%    effChans2Image={[1]}
%
%
%    %% option 2: provide all ids of channels 'nChannelsPerImage' & number
%    % of channel per effective channel 'nChannelsPerImage' channel
%    % example c: EXP 1: subcube 1 (first channel from each spectral
%    % window (2 of  which are flagged).
%    idChannels2Image = [1:16:272 289:16:320 337:16:512];
%    nChannelsPerImage = 1;
%
%    % example d: reduced imagecube containing 30 effective channels
%    % each concatenating 16 physical channels
%    idChannels2Image = [1:272 289:320 337:512];
%    nChannelsPerImage = 16;
%
% Warning
% -------
% - Data files: example of data filename: ``'data_ch_1.mat'``.
%
% - Any data file should contain the following variables : `maxProjBaseline`
%   (double, maximum projection basline), `u` , `v`, `w` (arrays of 
%   :math:`uvw`-coordinates), `nW` (noise whitening vector), `y` (vector of 
%   visibilities), `frequency` (associated frequency value).
%
% - Note that the data `y` are not whitened, uvw coordinates shoud be given
%   in units of the wavelength (i.e. normalised with the wavelength) and
%   `maxProjBaseline` in units of the wavelength.
%

% TODO give Matlab config file for a cluster (intead of just providing CIRRUS)
% TODO (keep empty parameter for non-used variables)
clear; clc; close all;
%% change paths if needed
% TODO: to be modified
main_dir = '..'; % '/Users/ad33/CodesScience/Faceted-Hyper-SARA/';
project_dir = [main_dir, filesep, 'experiments'];
cd(project_dir);

%% src name & datasets
% TODO: to be modified
imagecubeName = 'CYGA';
param_global.exp_type = 'test'; % ! only for debugging purposes
datasetNames = {'CYGA-ConfigA','CYGA-ConfigC'}; % allowing for multiple datasets
% data dir.
data_dir = [main_dir, filesep, 'data', filesep, imagecubeName, filesep];
fprintf('\nINFO: data are expected to be saved at %s\n', data_dir);
% preproc dir.
preproc_calib_dir = [data_dir, 'pre_processing_die/'];
% preproc_calib_dir=[data_dir,'pre_processing_dde/'];
% name of json parameter file
json_filename = "default_parameters.json";

% data files
% example of data file: 'data_ch_1.mat' with vars : 'maxProjBaseline', 'u','v','w', 'nW', 'y', 'frequency'.
% Note that data 'y' are not whitened, uvw coordinates are in units of the
% wavelength (i.e. normalised with the wavelength) and 'maxProjBaseline' (maximum projected baseline) is in
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
% > idChannels2Image = [1:16:272 289:16:320 337:16:512];
% > nChannelsPerImage = 1;

% example d: Exp 2: reduced imagecube containing 30 effective channels each concatenating 16 physical channels
% >idChannels2Image = [1:272 289:320 337:512];
% >nChannelsPerImage = 16;

% TODO: to be modified
idChannels2Image = 1:16; % ids of the 'physical' channels to be imaged
nChannelsPerImage = 16; % number of consecutive channels to be concatenated into each effective channel

nEffChans2Image = floor(numel(idChannels2Image) / nChannelsPerImage); % channels re-structured into effective channels
effChans2Image = cell(nEffChans2Image, 1);
for iEff = 1:nEffChans2Image
    if iEff < nEffChans2Image; effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:iEff * nChannelsPerImage);
    else; effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:end);
    end
    fprintf('\nEffective channel ID %d: physical channels involved: %d - %d\n', iEff, effChans2Image{iEff}(1), effChans2Image{iEff}(end));
end

%% running one subcube at a time
% TODO: to be modified
subcubeInd = 0; % id of subcube if spectral interleaving is active

% measurement op. params
param_global.measop_flag_dataReduction = 1;

% image details, dims &  cellsize
% TODO: to be modified
param_global.im_Nx = 2560; 
param_global.im_Ny = 1536; 
param_global.im_pixelSize =  0.06; % pixelsize in asec
param_global.algo_flag_computeOperatorNorm = true;

% faceting params: note that if interleaving is active, one subcube is
% imaged at a time: Qc=1 by default.
% TODO: to be modified
param_global.facet_Qx = 1; % dimFacet1
param_global.facet_Qy = 1; % dimFacet2
param_global.facet_overlap_fraction = [0.5, 0.5];

% reg params
% TODO: to be modified
param_global.reg_gam = 0.33; % l21 reg param
param_global.reg_gam_bar = 0.33; % nuclear norm reg param
param_global.reg_flag_reweighting = true;
param_global.reg_flag_homotopy = false; % not compulsory to be defined here..
param_global.reg_nReweights = 5;
param_global.generate_eps_nnls = false;

% algo & parallelisation params
% TODO: to be modified
param_global.algo_version = 'sara'; % 'fhs', 'hs' or 'sara'

% filenames and input
% TODO: to be modified
param_global.main_dir = main_dir;
% param_global.preproc_filename_dde = @(firstch,lastch) strcat(preproc_calib_dir,filesep,'ddes',filesep,'chs',num2str(firstch),'-',num2str(lastch),'_dies.mat');
param_global.preproc_filename_die = @(firstch, lastch) strcat(preproc_calib_dir, filesep,'dies',filesep,'chs', num2str(firstch), '-', num2str(lastch), '_dies.mat');
%param_global.preproc_filename_die = [];
%
param_global.preproc_filename_l2bounds = @(firstch, lastch) strcat(preproc_calib_dir, filesep,'l2bounds',filesep,'chs', num2str(firstch), '-', num2str(lastch), '_l2bounds.mat');
% param_global.preproc_filename_l2bounds = [];
%
param_global.preproc_filename_model = @(firstch, lastch) strcat(preproc_calib_dir, filesep,'model_images',filesep,'chs', num2str(firstch), '-', num2str(lastch), '_model_image.fits');
%param_global.preproc_filename_model = [];

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
