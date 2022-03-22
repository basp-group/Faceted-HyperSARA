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
%     containing pre-computed :math:`\ell_2` bounds. If not used, to be set
%     to ``[]``.
% param_global.preproc_filename_model : anonymous function
%     Function handle, taking the index of the first and last channel, and
%     returning the name of a file containing a model image to be used to
%     initialize the reconstruction algorithm. If not used, to be set to
%     ``[]``.
% param_global.measop_flag_visibility_gridding : bool
%     Flag to activate data dimensionality reduction via visbility gridding
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
%     Pixel-size in arcsec in the frequency domain. Set to ``[]`` to use
%     the default computation (implemented in
%     :mat:func:`imaging.imaging`).
% param_global.reg_flag_reweighting : bool
%     Flag to activate the use of multiple reweighting steps.
% param_global.reg_nReweights : int
%     Maximum number of reweighting steps.
% param_global.parcluster : string
%     String the name of the parcluster profile. By default ``"local"`` profile
%     is used. The user should set it to the name of the slurm parcluster
%     profile created on his/her HPC machine if prefered. 
% param_global.algo_solver : string
%     Name of the imaging approach used. Possible values are ``"fhs"``
%     (Faceted-HyperSARA), ``"hs"`` (HyperSARA) and ``"sara"`` (SARA).
% param_global.facet_Qx : int
%     Number of spatial facets along spatial axis x. Active only in
%     ``"fhs"``. If the field is not specified, the maximum possible
%     number of facets will the selected along the y axis.
% param_global.facet_Qy : int
%     Number of spatial facets along spatial axis y. Active only in
%     ``"fhs"``. If the field is not specified, the maximum possible
%     number of facets will the selected along the y axis.
% param_global.facet_overlap_fraction : array[double]
%     Array containing the overlap fraction between consecutive facets along
%     each axis (y and x) for the faceted low-rankness prior. Will be reset
%     automatically to ``[0, 0]`` if ``param_global.algo_solver = "sara"``
%     or ``"hs"``. Besides, each entry of
%     ``param_global.facet_overlap_fraction`` is reset to 0 the
%     number of facet along the corresponding dimension is equal to 1
%     (i.e., ``param_global.facet_Qy = 1`` and/or
%     ``param_global.facet_Qx = 1``).
% param_global.facet_subcubeInd : int
%     Index of the subcube to be reconstructed (if spectral faceting is
%     used). To be set to ``0`` otherwise.
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

%% Documentation
%
% For more details about the different parameters to be configured, see the
% online documentation: https://basp-group.github.io/Faceted-Hyper-SARA/default.html .

clear; clc; close all;

%% Parameters to be specified
% TODO give Matlab config file for a cluster (intead of just providing CIRRUS)
% TODO (keep empty parameter for non-used variables)

% change paths if needed
% TODO: to be adjusted by the user if required (keep current setting by default)
main_dir = '..';
project_dir = [main_dir, filesep, 'imaging'];
cd(project_dir);

% src name & datasets
% TODO: to be adjusted by the user
srcName = 'cygASband_Cube_256_512';
param_global.exp_type = 'test'; % ! only for debugging purposes
datasetNames = {'test'}; % allowing for multiple datasets, empty cell if a single dataset

% data directory
% TODO: to be adjusted by the user
data_dir = [main_dir, filesep, 'data', filesep, srcName, filesep];
fprintf('\nINFO: data are expected to be saved at %s\n', data_dir);
% preproc dir.
preproc_calib_dir = [data_dir, 'pre_processing_die/']; % set to [] if not used 
% preproc_calib_dir = [data_dir,'pre_processing_dde/'];
% name of json parameter file
json_filename = "default_parameters.json";

% data files
% example of data file: 'data_ch_1.mat' with vars : 'maxProjBaseline', 'u','v','w', 'nW', 'y', 'frequency'.
% Note that data 'y' are not whitened, uvw coordinates are in units of the
% wavelength (i.e. normalised with the wavelength) and 'maxProjBaseline' (maximum projected baseline) is in
% units of the wavelength
% ! Note the path for .mat files (dataset nametag used)
% TODO: to be adjusted by the user
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
% > idChannels2Image = [1:272 289:320 337:512];
% > nChannelsPerImage = 16;

% TODO: to be adjusted by the user
idChannels2Image = [1:2]; % ids of the 'physical' channels to be imaged, 
% e.g., [1:2] for channel range from 1 to 2
nChannelsPerImage = 1; % number of consecutive channels to be concatenated into each effective channel

nEffChans2Image = floor(numel(idChannels2Image) / nChannelsPerImage); % channels re-structured into effective channels
effChans2Image = cell(nEffChans2Image, 1);
for iEff = 1:nEffChans2Image
    if iEff < nEffChans2Image; effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:iEff * nChannelsPerImage);
    else; effChans2Image{iEff} = idChannels2Image((iEff - 1) * nChannelsPerImage + 1:end);
    end
    fprintf('\nEffective channel ID %d: physical channels involved: %d - %d\n', iEff, effChans2Image{iEff}(1), effChans2Image{iEff}(end));
end



% measurement op. params
% TODO: to be adjusted by the user
param_global.measop_flag_visibility_gridding = 0; % 1 if active, 0 otherwise: activating visibility gridding to reduce data dimensionality

% image details, dims &  cellsize
% TODO: to be adjusted by the user
param_global.im_Nx = 512; % 2560;
param_global.im_Ny = 256; % 1536;
param_global.im_pixelSize = []; % 0.06; % pixelsize in asec, use [] to use
% the default value set from the uv-coverage
param_global.algo_flag_computeOperatorNorm = true;

% faceting params: note that if interleaving is active, one subcube is
% imaged at a time: Qc=1 by default.
% TODO: to be adjusted by the user
param_global.facet_Qx = 1; % number of spatial facets along axis x (only used in Faceted HyperSARA)
param_global.facet_Qy = 1; % number of spatial facets along axis y (only used in Faceted HyperSARA)
param_global.facet_overlap_fraction = [0.5, 0.5];
param_global.facet_subcubeInd =0 ;% id of subcube if spectral interleaving is active, 0 if inactive

% reg params
% TODO: to be adjusted by the user
param_global.reg_gam = 1; % additional scaling factor for sparsity 
% regularization (using heuristics described in Thouvenin2021). Set 
% parameter to 1 by default.
param_global.reg_gam_bar = 1; % additional scaling factor for low-rankness regularization (using heuristics described in Thouvenin2021). Only active 
% for HyperSARA and Faceted HyperSARA. Set parameter to 1 by default.
param_global.reg_flag_reweighting = true;
param_global.reg_nReweights = 5;

% algo & parallelisation params
% TODO: to be adjusted by the user
param_global.algo_solver = 'fhs'; % 'fhs' (Faceted HyperSARA), 'hs' (HyperSARA) or 'sara' (SARA approach)

% filenames and input
% TODO: to be adjusted by the user (set path to [] whenever it is not used)
param_global.main_dir = main_dir;
param_global.preproc_filename_dde = [];
param_global.preproc_filename_die = [];
param_global.preproc_filename_l2bounds = [];
param_global.preproc_filename_model = [];
% param_global.preproc_filename_dde = @(firstch,lastch) strcat(preproc_calib_dir,filesep,'ddes',filesep,'chs',num2str(firstch),'-',num2str(lastch),'_dies.mat');
% param_global.preproc_filename_die = @(firstch, lastch) strcat(preproc_calib_dir, filesep,'dies',filesep,'chs', num2str(firstch), '-', num2str(lastch), '_dies.mat');
% param_global.preproc_filename_l2bounds = @(firstch, lastch) strcat(preproc_calib_dir, filesep,'l2bounds',filesep,'chs', num2str(firstch), '-', num2str(lastch), '_l2bounds.mat');
% param_global.preproc_filename_model = @(firstch, lastch) strcat(preproc_calib_dir, filesep,'model_images',filesep,'chs', num2str(firstch), '-', num2str(lastch), '_model_image.fits');


% Parcluster profile 
% TODO: to be adjusted by the user
param_global.parcluster = 'local'; % name of the parcluster profile to use.  If required, add 
% the name of  your own slurm parcluster profile 

%% read and set configuration from .json file
[param_global, param_solver, ...
    param_nufft, param_precond, dict] = ...
    read_json_configuration(json_filename, param_global);

%% run main job

% ---
% ! test: load ground truth image (for debugging purposes)
% exp_type = param_global.exp_type;
% switch exp_type
%     case "spatial"
%         image_name = 'cygASband_Cube_1024_2048_20';
%         spectral_downsampling = 1;
%         spatial_downsampling = 1;
%     case "spectral"
%         image_name = 'cygASband_Cube_256_512_100';
%         spectral_downsampling = 1;
%         spatial_downsampling = 1;
%     case "test"
%         image_name = 'cygASband_Cube_512_1024_20';
%         spectral_downsampling = 10;
%         spatial_downsampling = 2;
%     case "local_test"
%         image_name = 'cygASband_Cube_256_512_100';
%         spectral_downsampling = 25;
%         spatial_downsampling = 1;
%         coverage_path = "data/vla_7.95h_dt10s.uvw256.mat";
%     case "old_local_test"
%         image_name = 'cubeW28';
%         spectral_downsampling = 20;
%         spatial_downsampling = 4;
%         coverage_path = "data/vla_7.95h_dt10s.uvw256.mat";
%     otherwise
%         error("Unknown experiment type");
% end
% reference_cube_path = fullfile('../data', strcat(image_name, '.fits'));
% info        = fitsinfo(reference_cube_path);
% rowend      = info.PrimaryData.Size(1);
% colend      = info.PrimaryData.Size(2);
% sliceend    = info.PrimaryData.Size(3);
% if strcmp(param_global.algo_solver, 'sara')
%     x0 = fitsread(reference_cube_path, 'primary', ...
%         'Info', info, ...
%         'PixelRegion', {[1 spatial_downsampling rowend], ...
%         [1 spatial_downsampling colend], ...
%         ind});
% else
%     x0 = fitsread(reference_cube_path, 'primary', ...
%         'Info', info, ...
%         'PixelRegion', {[1 spatial_downsampling rowend], ...
%         [1 spatial_downsampling colend], ...
%         [1 spectral_downsampling sliceend]});
% end
% [Ny, Nx, nchans] = size(x0);
% N = Nx * Ny;
% x0 = reshape(x0, [N, nchans]);
% ---

imaging(srcName, datasetNames, dataFilename, ...
     effChans2Image, param_solver, param_global, ...
    param_nufft, param_precond, dict);
