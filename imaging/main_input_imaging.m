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
%     data extraction, used to get the data and generate output files.
% datasetNames : cell of string
%     Name of the different datasets used (e.g., in the cases of
%     different acquisition configurations or different observation times),
%     to be set to ``[]`` if one data set is imaged & no name tag of
%     the MS is given during data extraction.
% dataFilename : anonymous function
%     Function handle taking a channel index as an input, and returning the
%     name of the associated data ``.mat`` file (one file per frequency channel
%     and dataset) following the nomenclature adopted during data extraction.
%     Do not modify unless the nomenclature has been altered.
% idChannels2Image : array[int]
%     Indices of the  input (`physical`) channels to be imaged.
% nChannelsPerImage : int
%     Number of consecutive channels to be concatenated into each  output
%     (effective) channel (set to `1` if full spectral resolution is
%     considered, i.e. number of output channels is equal to the number of
%     input channels).
% effChans2Image: cell array
%     Cell array containing the ids of the input (physical) channels to be
%     concatenated for each output (effective) channel. Set automatically if
%     `idChannels2Image` and `nChannelsPerImage` are provided.
% main_dir : string
%     Directory of the Faceted-HyperSARA repository.
% data_dir : string
%     Directory of the input ``.mat`` data files. Default ``./Faceted-HyperSARA/data``.
% imaging_dir : string
%     Directory of the imaging experiment. Default
%     ``./Faceted-HyperSARA/imaging``.
% calib_type: string
%     calibration type: `dde`, `die` or `[]` if no results from a calib
%     pre-processing step.
% preproc_calib_dir : string
%     Sub-directory in `data_dir`  containing files from a calibration
%     pre-processing step. To be set to ``[]`` if not used or not
%     available.
% json_filename : string
%     Name of the input ``.json`` configuration file specifying the
%     value of the algorithm and measurement operators' parameters
%     (PDFB, reweighting, ellipsoid,
%     projection, preconditioning, w-projection, noise estimation).
% param_global.im_Nx  : int
%     Number of pixels along axis x in the reconstructed image.
% param_global.im_Ny : int
%     Number of pixels along axis y in the reconstructed image.
% param_global.im_pixelSize  : double
%     Pixel-size in arcsec. Set to ``[]`` to use
%     the default value corresponding to 2 times the resolution of the
%     observation (given by the longest baseline at the highest frequency).
% param_global.adjust_flag_noise : bool
%     Flag to activate noise estimation on the fly in the case of
%     unreliable noise statistics.
% param_global.measop_flag_visibility_gridding : bool
%     Flag to activate visibility gridding for data dimensionality
%     reduction.
% param_global.data_flag_apply_imaging_weights : bool
%     Flag to read and apply imaging weights (e.g. Briggs, uniform) to
%     data. Typically, these are stored in the ``MS`` as ``IMAGING_WEIGHTS`` or
%     ``IMAGING_WEIGHTS_SPECTRUM``. If one of these columns
%     are found, imaging weights will be automatically extracted as ``nWimag``
%     in the ``.mat`` data file.
% param_global.algo_solver : string
%     Name of the imaging approach used. Possible values are ``"fhs"``
%     (Faceted-HyperSARA), ``"hs"`` (HyperSARA) and ``"sara"`` (SARA).
%     Default: monochromatic imaging ``"sara"``, wideband imaging ``"fhs"``.
% param_global.facet_Qx : int
%     Number of spatial facets along spatial axis x. Active only in
%     ``"fhs"``. If the field is not specified, the maximum possible
%     number of facets will the selected along the y axis.
% param_global.facet_Qy : int
%     Number of spatial facets along spatial axis y. Active only in
%     ``"fhs"``. If the field is not specified, the maximum possible
%     number of facets will the selected along the y axis.
% param_global.facet_overlap_fraction : double[2]
%     Fraction of the total size of facet overlapping with a neighbouring
%     facet along each axis (y and x) for the faceted low-rankness prior.
%     Active only in ``"fhs"``. Will be reset
%     automatically to ``[0, 0]`` if the spatial faceting is not active in
%     ``"fhs"`` along at least one dimension
%     (i.e., ``param_global.facet_Qy = 1`` and/or ``param_global.facet_Qx = 1``).
% param_global.facet_subcubeInd : int
%     Index of the subcube to be reconstructed (if spectral faceting is
%     used). To be set to ``0`` otherwise.
% param_global.reg_gam : double
%     Additional multiplicative factor affecting the joint average sparsity
%     regularization term in ``"fhs"`` and ``"hs"``. Additional
%     multiplicative factor affecting the average sparsity
%     regularization term in ``"sara"``.
% param_global.reg_gam_bar : double
%    Additional multiplicative factor affecting  low-rankness prior
%    regularization parameter. Active only in ``"fhs"`` and ``"hs"``.
% param_global.reg_flag_reweighting : bool
%     Flag to activate the re-weighting procedure. Default set to
%     ``true``.
% param_global.reg_nReweights : int
%     Maximum number of reweighting iterations.
% param_global.parcluster : string
%     Name of the parallel parcluster profile to launch the parpool. By
%     default  ``"local"`` profile is used. The user should set it to the
%     name of the slurm parcluster profile created on his/her HPC machine,
%     if prefered.
% param_global.algo_flag_computeOperatorNorm (bool)
%     Flag to activate computation of the measurement operator norm,
%     default, true.
% param_global.algo_flag_saveOperatorNorm (bool)
%     Flag to activate saving of the measurement operator norm,
%     default, false.
% param_global.algo_flag_solveMinimization (bool)
%     Flag to trigger the solver, true.
% param_global.preproc_filename_noise_std : anonymous function
%     Function handle, taking the indices of the physical channel
%     and the dataset, and returning a string corresponding to the name of a file
%     containing noise statistics obtained from a pre-processing step.
%     Expected vars: ``sigma``: the standard deviation of the noise,
%     ``RESIDUAL``: (optional) the residual visibilities from a pre-processing
%     step. If not used, should be set to ``[]``.
% param_global.preproc_filename_cal_solutions : anonymous function
%     Function handle, taking the indices of the physical channel and the
%     dataset, and returning a string corresponding to the name of a file
%     containing DDE/DIE calibration solutions.
%     Expected var: ``DDEs``: complex array if DDE calibration and
%     ``DIEs``: complex vector if DIE calibration . If not used, should be set to ``[]``.
% param_global.preproc_filename_model : anonymous function
%     Function handle, taking the indices of the first and last physical channels,
%     associated with an effective (output) channel, and
%     returning the name of a file containing a model image to be used to
%     initialize the reconstruction algorithm. If not used, should be set to
%     ``[]``.
%
% Example
% -------
%
% .. code-block:: matlab
%
%    %% Examples to set the input (physical) and output (effective) channels IDs
%    %% Option 1: provide  `effChans2Image` a cell array containing the ids
%    % of the physical channels to be concatenated for each effective channel.
%    % example a: two effective channels, containing two physical
%    % channels each
%    effChans2Image = {[1,2], [3,4]};
%
%    % example b: one effective channel with one physical channel
%    effChans2Image = {[1]};
%
%
%    %% Option 2: provide all ids of physical channels `nChannelsPerImage` & number
%    % of channel per effective channel `nChannelsPerImage`. The list of the ids
%    % of the physical channels to be concatenated for each effective
%    % channel `effChans2Image` will be created automatically.
%    % example c:  subcube 1 (first channel from each spectral
%    % window (in this example two spectral windows are flagged).
%    % The number of output (effective) channels is equal to the number
%    % of input (physical) channels
%    idChannels2Image = [1:16:272, 289:16:304, 321:16:512];
%    nChannelsPerImage = 1;
%
%    % example d: reduced imagecube containing 30 effective channels
%    % each combining 16 physical channels.
%    % In this example the number of output (effective) channels is reduced to 30 from
%    % 480 input (physical) channels.
%    idChannels2Image = [1:272, 289:304, 321:512];
%    nChannelsPerImage = 16;
%
% Warning
% -------
% - Data files (example: ``'data_ch_1.mat'``) are expected to contain the following variables : ``y`` (complex vector of
%   visibilities, Stokes I), ``u`` , ``v``, ``w`` (arrays of
%   :math:`uvw`-coordinates), `maxProjBaseline`
%   (double, maximum projected baseline), ``nW`` (noise whitening vector),
%   ``frequency`` (associated frequency) and ``nWimag`` (optional,
%   imaging weights such as uniform of Briggs)..
% - Extracted data ``y``  are not whitened (i.e. natural weights not applied).
% - The variables ``maxProjBaseline`` and the :math:`uvw` coordinates are
%   in units of the wavelength (i.e. normalised by the associated wavelength).
% - When spectral faceting is enabled, run ``imaging/main_input_imaging.m`` for each
%   sub-cube, by updating the ID of the
%   sub-cube to image ``param_global.facet_subcubeInd`` along with the
%   associated channel frequencies ``effChans2Image``.
%

%% Documentation
%
% For more details about the different parameters to be configured, see the
% online documentation: https://basp-group.github.io/Faceted-HyperSARA/default.html .

clear all; close all;
%% ------------------------------------------------------------------------
%% Parameters to be specified
% / TODO (keep empty parameter for non-used variables)

%% target source & datasets names as specified during data extraction
% / TODO: to be adjusted by the user
srcName = 'CYGA'; % src tag as specified during data extraction, {''} otherwise.
datasetNames = {'CYGA-ConfigA', 'CYGA-ConfigC'}; % accomodating multiple datasets,
% set as specified in the data extraction step, [] otherwise.

%% Paths
% Main path of the Faceted-HyperSARA repository
% / TODO: to be adjusted by the user if required.
% Keep current setting by default
main_dir = ['..', filesep];

% Imaging experiment directory
imaging_dir = [main_dir, filesep, 'imaging'];
cd(imaging_dir);

% data directory, as expected from the data extraction step
data_dir = [main_dir, filesep, 'data', filesep, srcName, filesep];
fprintf('\nINFO: data are expected to be saved at %s\n', data_dir);

% data files
% ! To be adjusted only if data paths during data extraction have been altered
% example of a data file: `./data/CYGA/CYGA-ConfigA/data_ch_1.mat`
% ! Note in the path for .mat files, dataset name tag are used.
dataFilename = @(idSet, ch) strcat(data_dir, datasetNames{idSet}, filesep, 'data_ch_', num2str(ch), '.mat');

% calibration preprocessing directory
% ! set path to [] whenever it is not used
preproc_calib_dir = [data_dir, 'pre_processing', filesep];
preproc_results_exist = isfolder(preproc_calib_dir);

%% read and set configuration from .json file
json_filename = "default_parameters.json"; % json parameter file
[param_global, param_solver, param_nufft, param_precond, param_wproj, dict] = ...
    read_json_configuration(json_filename);
    
%% Faceted-HyperSARA path
param_global.main_dir = main_dir;

%% image details, dims &  cellsize
% / TODO: to be adjusted by the user
param_global.im_Nx = 2560;
param_global.im_Ny = 2560;
param_global.im_pixelSize =  0.06; % pixelsize in asec, set to [] to use the default value set from the uv-coverage.

%% channel organisation
% !! option 1: provide a cell array containing the ids of the  channels to be concatenated for each effective channel.
% -- example a: two effective (output) channels, combining two 'physical'
% (input) channels each
% > effChans2Image={[1,2],[3,4]};
% -- example b: one effective channel with one physical channel
% > effChans2Image={[1]};
% !! option 2: provide all ids of channels to be imaged 'nChannelsPerImage' & num of channel per effective channel 'nChannelsPerImage'.
% 'effChans2Image' will then be determined automatically.
% -- example c: EXP 1: subcube 1 (first channel from each spectral window (2 of  which are flagged).
% > idChannels2Image = [1:16:272 289:16:320 337:16:512];
% > nChannelsPerImage = 1;
% -- example d: Exp 2: reduced imagecube containing 30 effective channels each concatenating 16 physical channels
% > idChannels2Image = [1:272 289:320 337:512];
% > nChannelsPerImage = 16;
% / TODO: to be adjusted by the user
idChannels2Image =  1:256; % ids of the physical (input) channels to be imaged, e.g., [1:256] for channel range from 1 to 256
nChannelsPerImage = 16; % number of consecutive channels to be combined into each effective (output) channel.

% re-arranging channels indices based on idChannels2Image, nChannelsPerImage
nEffChans2Image = floor(numel(idChannels2Image) / nChannelsPerImage); % num of ouput effective channels: number of images in the estimate model cube
effChans2Image = cell(nEffChans2Image, 1);
for iEff = 1:nEffChans2Image
    if iEff < nEffChans2Image
        effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:iEff * nChannelsPerImage);
    else
        effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:end);
    end
    fprintf('\nINFO:Effective channel ID %d: physical channels involved: %d - %d', ...
        iEff, effChans2Image{iEff}(1), effChans2Image{iEff}(end));
end

%% measurement operator features
% / TODO: to be adjusted by the user
% activating visibility gridding for data dimensionality reduction
param_global.measop_flag_visibility_gridding = true; % 1 if active, 0 otherwise.

%% noise level estimation
% / TODO: to be adjusted by the user: should be set to true if reliable noise estimates are not available
param_global.adjust_flag_noise = true;

%% algo params
% / TODO: to be adjusted by the user:
% 'fhs' (Faceted HyperSARA), 'hs' (HyperSARA) or 'sara' (SARA approach)
param_global.algo_solver = 'fhs';

%% faceting params
%  note that if spectral faceting is active, one subcube (i.e. spectral facet) is imaged at a time Qc=1 by default.
% / TODO: to be adjusted by the user
param_global.facet_Qx = 4; % num. of spatial facets along axis x (only used in Faceted HyperSARA)
param_global.facet_Qy = 4; % num. of spatial facets along axis y (only used in Faceted HyperSARA)
param_global.facet_overlap_fraction = [0.1, 0.1]; % (only used in Faceted HyperSARA)
param_global.facet_subcubeInd = 0; % id of the subcube to image if spectral faceting is active, 0 otherwise

%% regularisation params
% / TODO: to be adjusted by the user
param_global.reg_gam = 1; % default 1, additional scaling factor for
% sparsity regularization (using heuristics described in Thouvenin2021).
param_global.reg_gam_bar = 1; % default 1, additional scaling factor for
% low-rankness regularization (using heuristics described in Thouvenin2021).
% Active only for HyperSARA and Faceted HyperSARA.
param_global.reg_flag_reweighting = true;
param_global.reg_nReweights = 30;

%% filenames and input from a pre-processing step
if preproc_results_exist
    % / TODO: to be adjusted by the user
    % ! set path to [] whenever it is not used
    % noise level files
    param_global.preproc_filename_noise_std = @(idSet, ch) strcat(preproc_calib_dir, ...
         'sigma',filesep,datasetNames{idSet}, filesep, 'chs', num2str(ch), '_sigma.mat');

     % init. images files
    param_global.preproc_filename_model =  @(firstch, lastch) strcat(...
         preproc_calib_dir, 'model_images', filesep, 'chs', num2str(firstch),...
         '-', num2str(lastch), '_model_image.fits');

    % ! dde calib / die calib files
    param_global.preproc_filename_cal_solutions = @(idSet, ch) strcat(preproc_calib_dir, ...
        'cal_solutions', filesep, datasetNames{idSet}, filesep, 'chs', num2str(ch), ...
        '_cal_solutions.mat');
end

%% Parallel parcluster profile
% / TODO: to be adjusted by the user
% name of the parallel parcluster profile to use. To be adjusted if running on a HPC
% and a slurm parcluster profile will be used. Default set to 'local'.
param_global.parcluster = 'mySlurmProfileSingleThread';

%% ------------------------------------------------------------------------
%% run  imaging job
imaging(srcName, datasetNames, dataFilename, effChans2Image, param_solver, ...
    param_global,  param_nufft, param_precond, param_wproj, dict);
