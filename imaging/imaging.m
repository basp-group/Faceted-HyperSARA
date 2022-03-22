function imaging(srcName, datasetsNames, dataFilename, ...
     effChans2Image, param_solver, param_global, ...
    param_nufft, param_precond, dict)
% Main function to run any algorithm from the SARA family using the
% configuration provided in :mat:scpt:`imaging.main_input_imaging`.
% The approaches are Faceted HyperSARA, HyperSARA and SARA.
%
% Parameters
% ----------
% srcName : string
%     Name of the target source to be reconstructed as specified during
%     data extraction, used to get the data and generate output files.
% dataSetsNames: cell of string
%     Names of the different datasets used (e.g., in the cases of
%     different acquisition configurations or different observation times),
%     to be set to ``{''}`` if one data set is imaged or no name tag of 
%     the MS is given during data extraction. 
% dataFilenames: cell of string
%     Function handle taking a channel index as an input and returning the
%     name of the associated data ``.mat`` file (one file per frequency channel),
%     following the nomenclature adopted during data extraction .
% effChans2Image: cell array
%     Cell array containing the ids of the physical (input) channels to be
%     concatenated for each effective (output) channel.
% param_global: struct
%     Global imaging pipeline parameters (see
%     :mat:func:`imaging.main_input_imaging`).
% param_global.im_pixelSize  : double
%     Pixel-size in arcsec. Set to ``[]`` to use
%     the default value corresponding to 2 times the resolution of the
%     observation (given by the longest baseline at the highest frequency).
% param_global.im_Nx : int
%     Image dimension (x axis, dimension 2).
% param_global.im_Ny : int
%     Image dimension (y axis, dimension 1).
% param_global.measop_flag_visibility_gridding : bool
%     Flag to activate data dimensionality reduction via visibility gridding in the
%     definition of the measurement operator.
% param_global.algo_solver : string (``"sara"``, ``"hs"`` or ``"fhs"``).
%     Selected solver.
% param_global.facet_subcubeInd: int
%     Index of the spectral facet to be reconstructed (set to -1 or 0 to
%     deactivate spectral faceting).
% param_global.facet_Qc : int
%     Number of spectral facets. Active only in ``"fhs"``.
% param_global.facet_Qx : int
%     Number of spatial facets along spatial axis x. Active only in
%     ``"fhs"``.
% param_global.facet_Qy : int
%     Number of spatial facets along spatial axis y. Active only in
%     ``"fhs"``.
% param_global.facet_overlap_fraction : double[2]
%     Fraction of the total size of facet overlapping with a neighbouring
%     facet along each axis (y and x) for the faceted low-rankness prior. 
%     Active only in ``"fhs"``. Will be reset
%     automatically to ``[0, 0]`` if the spatial faceting is not active in ``"fhs"``
%     along at least one dimension(i.e., ``param_global.facet_Qy = 1`` and/or ``param_global.facet_Qx = 1``).
% param_global.facet_window_type : string (``"triangular"``, ``"hamming"`` or ``"pc"`` (piecewise-constant))
%     Type of apodization window considered for the faceted nuclear norm
%     prior. Active only in ``"fhs"``.
% param_global.reg_gam : double
%     Additional multiplicative factor affecting the joint average sparsity
%     regularization term in ``"fhs"`` and ``"hs"``. Additional 
%     multiplicative factor affecting the average sparsity
%     regularization term in ``"sara"``.
% param_global.reg_gam_bar : double
%    Additional multiplicative factor affecting  low-rankness prior 
%    regularization parameter. Active only in ``"fhs"`` and ``"hs"``.
% param_global.reg_flag_reweight : int
%     Flag to activate re-weighting.
% param_global.reg_nReweights : int
%     Maximum number of reweighting steps.
% param_global.algo_flag_computeOperatorNorm : bool
%     Flag to trigger the computation of the norm of the (preconditionned) measurement
%     operator. Default set to ``true``. If set to ``false``,
%     MATLAB will look for a file where this
%     quantity has been saved (save and computation steps are triggered in
%     :mat:func:`imaging.imaging`).
% param_global.algo_flag_solveMinimization : bool
%     Flag triggering the solver (``"fhs"``, ``"hs"`` or ``"sara"``).
% param_global.ncores_data : int
%     Number of cores handlig the data fidelity terms (data cores). For
%     Faceted HyperSARA, the total number of cores used
%     is ``Qx*Qy + ncores_data + 1``. For SARA and HyperSARA, represents
%     the number of cores used for the parallelization.
% param_global.parcluster : string
%     Name of the parallel parcluster profile to launch the parpool. By default ``"local"`` profile
%     is used. The user should set it to the name of the slurm parcluster
%     profile created on his/her HPC machine, if prefered. 
% param_global.preproc_filename_l2bounds : anonymous function
%     Function handle, taking the index of the first and last physical channels, and
%     returning a string corresponding to the name of a file
%     containing pre-computed :math:`\ell_2` bounds. If not used, must be set
%     to ``[]``.
% param_global.preproc_filename_model : anonymous function
%     Function handle, taking the index of the first and last physical channels, and
%     returning the name of a file containing a model image to be used to
%     initialize the reconstruction algorithm. If not used, should be set to
%     ``[]``.
% param_global.preproc_filename_die : anonymous function
%     Function handle, taking the index of the first and last physical channels, and
%     returning a string corresponding to the name of a file containing
%     DIE calibration constants. If not used, must be set to ``[]``.
% param_global.preproc_filename_dde : anonymous function
%     Function handle, taking the index of the first and last physical channels, and
%     returning a string corresponding to the name of a file containing
%     DDE calibration constants. If not used, must be set to ``[]``.

% Note
% ----
% - Examples to set the input (physical) and output (effective) channels IDs:
%
% .. code-block:: matlab
%
%     %% Name of the dataset to be imaged
%     % example a: one dataset -name tag is not compulsory:
%     datasetNames = {''};
%     % example b: two data sets from two configurations of the VLA
%     datasetNames = {'CYGA-ConfigA','CYGA-ConfigC'};
%
%     %% Indices of the 'physical' channels to be concatenated
%     % example a: two effective channels, containing two 'physical'
%     % channels each
%     effChans2Image = {[1,2], [3,4]};
%     % example b: one channel effective channel with one physical channel
%     effChans2Image = {[1]};
%
%

%
% Deprecated fields
%
% param_global.reg_flag_homotopy : int
%     Flag to activate homotopy strategy in the re-weighting scheme.
%

% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
%% constant
speed_of_light = 299792458;
facet_dim_min  = 256;
default_superresolution = 2;
%% check input params : @PA in light of the new structure of params, many of these fields will be deleted I guess ..
% Image resolution & dimensions
% ! parameters and default values (keep only what is required)
if ~isfield(param_global, 'im_Nx');  param_global.im_Nx = 2048; end
if ~isfield(param_global, 'im_Ny');  param_global.im_Ny = 2048; end
if ~isfield(param_global, 'im_pixelSize');  param_global.im_pixelSize = []; end

% Prior: wideband facet related
if ~isfield(param_global, 'facet_subcubeInd');  param_global.facet_subcubeInd = 0; end
if ~isfield(param_global, 'facet_Qx');  param_global.facet_Qx =  floor(param_global.im_Nx / facet_dim_min);
else
    facet_dim = floor(param_global.im_Nx / param_global.facet_Qx );
    if facet_dim < facet_dim_min
        fprintf('\nWARNING: facet dimension is small, advised max. value:  param_global.facet_Qx=%d.', floor(param_global.im_Nx / facet_dim_min));
    end
end
if ~isfield(param_global, 'facet_Qy');  param_global.facet_Qy =  floor(param_global.im_Ny / facet_dim_min);
    facet_dim = floor(param_global.im_Ny / param_global.facet_Qy );
    if facet_dim < facet_dim_min
        fprintf('\nWARNING: facet dimension is small, advised max. value: param_global.facet_Qy=%d.', floor(param_global.im_Ny / facet_dim_min));
    end
end
if ~isfield(param_global, 'facet_Qc'); param_global.facet_Qc = 1; end
if ~isfield(param_global, 'facet_window_type'); param_global.facet_window_type = 'triangular'; end
if ~isfield(param_global, 'facet_overlap_fraction'); param_global.facet_overlap_fraction = [0.2, 0.2]; end


% Prior: sparsity dict.
if ~isfield(dict, 'nlevel'); dict.nlevel = 4;  end
if ~isfield(dict, 'basis'); dict.basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};  end
% Prior: reg params
if ~isfield(param_global, 'reg_gam'); param_global.reg_gam = 1;  end
if ~isfield(param_global, 'reg_gam_bar'); param_global.reg_gam_bar = 1;  end
if ~isfield(param_global, 'reg_flag_reweighting'); param_global.reg_flag_reweighting = 1; end


if param_global.reg_flag_reweighting &&  ~isfield(param_global, 'reg_nReweights')
    param_global.reg_nReweights = 5;
end
% if ~isfield(param_global, 'reg_flag_homotopy'); param_global.reg_flag_homotopy = 0; end
% Data blocks
% if ~isfield(param_global, 'generate_eps_nnls'); param_global.generate_eps_nnls = false; end
% if ~isfield(param_global, 'data_nDataBlk');  param_global.data_nDataBlk = []; end
% if ~isfield(param_global, 'data_sizeDataBlk');  param_global.data_sizeDataBlk = []; end
% if isempty(param_global.data_sizeDataBlk) && isempty (param_global.data_nDataBlk)
%     param_global.data_sizeDataBlk = 2e5;
% end

% Algo
if ~isfield(param_global, 'algo_solver')
    if numel(effChans2Image)>1,  param_global.algo_solver = 'fhs';
    else, param_global.algo_solver = 'sara';
    end
end
if ~isfield(param_global, 'algo_flag_solveMinimization'); param_global.algo_flag_solveMinimization = true; end
if ~isfield(param_global, 'algo_flag_computeOperatorNorm'); param_global.algo_flag_computeOperatorNorm = true; end
if ~isfield(param_global, 'algo_ncores_data'); param_global.algo_ncores_data = numel(effChans2Image);  end

% Measurement operator
if ~isfield(param_global, 'measop_flag_visibility_gridding'); param_global.measop_flag_visibility_gridding = false; end
if ~isfield(param_global, 'measop_flag_wproj'); param_global.measop_flag_wproj = []; end
if param_global.measop_flag_wproj
    if ~isfield(param_global, 'measop_wprojCEnergyL2'); param_global.measop_wprojCEnergyL2 = 1 - 1e-4; end
    if ~isfield(param_global, 'measop_wprojGEnergyL2'); param_global.measop_wprojGEnergyL2 = 1 - 1e-4; end
else
    param_global.measop_wprojCEnergyL2 = 1;
    param_global.measop_wprojGEnergyL2 = 1;
end

% Project dir.
if ~isfield(param_global, 'main_dir'); param_global.main_dir = [pwd, filesep]; end

% Input filenames
if ~isfield(param_global, 'preproc_filename_model'); param_global.preproc_filename_model = []; end
if ~isfield(param_global, 'preproc_filename_l2bounds'); param_global.preproc_filename_l2bounds = []; end
if ~isfield(param_global, 'preproc_filename_die'); param_global.preproc_filename_die = []; end
if ~isfield(param_global, 'preproc_filename_dde'); param_global.preproc_filename_dde = []; end
% Output filenames
if ~isfield(param_global, 'exp_type'); param_global.exp_type = '';  end %  AD: maybe remove?

% name of the parallel parcluster profile
if ~isfield(param_global, 'parcluster'); param_global.parcluster = 'local';
end 


% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% get params
% image dimensions & resolution
Nx = param_global.im_Nx;
Ny = param_global.im_Ny;
pixelSize  = param_global.im_pixelSize;
if isempty (pixelSize)
    imResolution = 'nominal';
else
    imResolution = 'user_defined';
end
switch imResolution
    case 'nominal'
        fprintf('\nWARNING: pixelsize not found, default is considered.\n');
    otherwise
        fprintf('\nINFO: pixelsize: %f arcsec.\n', pixelSize);
end
% Prior: wideband facet related
subcubeInd = param_global.facet_subcubeInd;
Qx = param_global.facet_Qx;
Qy = param_global.facet_Qy;
Qc = param_global.facet_Qc;
window_type = param_global.facet_window_type;
overlap_fraction = param_global.facet_overlap_fraction;
% Prior: reg params
gam = param_global.reg_gam;
gam_bar =  param_global.reg_gam_bar;
flag_reweighting = param_global.reg_flag_reweighting;
if flag_reweighting
    reg_nReweights = param_global.reg_nReweights;
end
% Measurement operator settings
flag_visibility_gridding = param_global.measop_flag_visibility_gridding; % data reduction
param_wproj.do = param_global.measop_flag_wproj; % w projection
if isempty(param_global.measop_flag_wproj)
    param_wproj.do = 0;
end
param_wproj.CEnergyL2 = param_global.measop_wprojCEnergyL2; % w projection
param_wproj.GEnergyL2 = param_global.measop_wprojGEnergyL2; % w projection
% Algo
algo_solver = param_global.algo_solver;
ncores_data = param_global.algo_ncores_data;
flag_computeOperatorNorm = param_global.algo_flag_computeOperatorNorm;
flag_solveMinimization =  param_global.algo_flag_solveMinimization;
% Pre-processing step
param_preproc.filename_die  = param_global.preproc_filename_die;
param_preproc.filename_dde  = param_global.preproc_filename_dde;
param_preproc.filename_l2bounds = param_global.preproc_filename_l2bounds;
param_preproc.filename_model = param_global.preproc_filename_model;
param_preproc.subcube = subcubeInd;
param_preproc.done = (~isempty(param_preproc.filename_l2bounds)) * ~isempty(param_preproc.filename_model) * (~isempty(param_preproc.filename_die) || ~isempty(param_preproc.filename_dde));
filename_l2bounds = param_global.preproc_filename_l2bounds; % available l2bounds
flag_calib.die = 0;
flag_calib.dde = 0;
if ~isempty(param_preproc.filename_die)
    flag_calib.die = 1;
elseif ~isempty(param_preproc.filename_dde)
    flag_calib.dde = 1;
end

% parcluster
parcluster = param_global.parcluster;
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% Paths
format compact;
project_dir = param_global.main_dir;
fprintf('\nMain project dir. is %s: ', project_dir);
current_dir = pwd;
if strcmp(project_dir, current_dir)
    current_dir = [project_dir, 'imaging', filesep];
    cd(current_dir);
end
fprintf('\nCurrent dir. is  %s: ', current_dir);
% src & lib codes
addpath([project_dir, filesep, 'lib', filesep, 'operators', filesep]);
addpath([project_dir, filesep, 'lib', filesep, 'utils', filesep]);
addpath([project_dir, filesep, 'lib', filesep, 'measurement-operator', filesep, 'nufft']);
addpath([project_dir, filesep, 'lib', filesep, 'measurement-operator', filesep, 'lib', filesep, 'utils']);
addpath([project_dir, filesep, 'lib', filesep, 'measurement-operator', filesep, 'lib', filesep, 'operators']);
addpath([project_dir, filesep, 'lib', filesep, 'measurement-operator', filesep, 'lib', filesep, 'ddes_utils']);
addpath([project_dir, filesep, 'lib', filesep, 'faceted-wavelet-transform', filesep, 'src']);
addpath([project_dir, filesep, 'src']);
addpath([project_dir, filesep, 'src', filesep, 'heuristics', filesep]);
%
if strcmp(algo_solver, 'sara')
    addpath([project_dir, filesep, 'src', filesep, 'sara']);
elseif strcmp(algo_solver, 'hs')
    addpath([project_dir, filesep, 'src', filesep, 'hs', filesep]);
elseif strcmp(algo_solver, 'fhs')
    addpath([project_dir, filesep, 'src', filesep, 'fhs', filesep]);
end
% setting paths to results and reference image cube
results_path = fullfile(project_dir,'results', srcName);
mkdir(results_path);
auxiliary_path = fullfile(results_path, algo_solver);
mkdir(auxiliary_path);
fprintf('\nINFO: results will be saved in: \n %s .\n', auxiliary_path);
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% info
disp('Imaging configuration');
disp(['Source name: ', srcName]);
disp(['Image cube size: ', num2str(Ny), ' x ', num2str(Nx), ' x ', num2str(numel(effChans2Image))]);
disp(['Algorithm version: ', algo_solver]);
disp(['Number of facets Qy x Qx : ', num2str(Qy), ' x ', num2str(Qx)]);
disp(['Number of spectral facets Qc : ', num2str(Qc)]);
if ~strcmp(algo_solver, 'sara')
    disp(['Overlap fraction: ', strjoin(strsplit(num2str(overlap_fraction)), ', ')]);
end
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% global params
N = Nx * Ny;
nEffectiveChans  = numel(effChans2Image);
ImageCubeDims = [Ny, Nx, nEffectiveChans];
nDataSets = numel(datasetsNames);
% set pixelsize
if strcmp(imResolution, 'nominal')
    maxProjBaseline = 0;
    for iEffCh = 1:nEffectiveChans
        for iCh = 1:numel(effChans2Image{iEffCh})
            for idSet = 1:nDataSets
                dataloaded = load(dataFilename(idSet, effChans2Image{iEffCh}(iCh)), 'maxProjBaseline');
                maxProjBaseline = max(maxProjBaseline, dataloaded.maxProjBaseline);
            end
        end
    end; clear dataloaded;
    spatialBandwidth = 2 * maxProjBaseline;
    pixelSize = (180 / pi) * 3600 / (default_superresolution * spatialBandwidth);
    fprintf('\nINFO: default pixelsize: %f arcsec, that is %d x nominal resolution at the highest freq.', default_superresolution,pixelSize);
end
halfSpatialBandwidth = (180 / pi) * 3600 / (pixelSize) / 2;
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% image cube initialisation: load if available otherwise set to 0
switch algo_solver
    case 'sara'
        xinit = zeros(Ny, Nx);
        if ~isempty(param_preproc.filename_model)
            fprintf('\nLoading  available image estimate  ...');
            try xinit = fitsread(param_preproc.filename_model(effChans2Image{1}(1), effChans2Image{1}(end)));
                if  ~(size(xinit, 1) == Ny && size(xinit, 2) == Nx)
                    xinit = zeros(Ny, Nx);
                    fprintf('\nWARNING: init. image not found.');
                else; fprintf('\nINFO: init. image found.');
                end
            catch ; fprintf('\nWARNING: init. image not found.');
            end
        end
    otherwise
        xinit = zeros(N, nEffectiveChans);
        if ~isempty(param_preproc.filename_model)
            fprintf('\nLoading  available image estimate  ...');
            for  iEffCh =  1:nEffectiveChans
                try  xinit(:, iEffCh) = reshape(fitsread(param_preproc.filename_model(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end))), N, 1);
                    fprintf('\nINFO: init. image found.');
                catch ; fprintf('\nWARNING: init. image not found.');
                end
            end
        end
end
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% l2 bounds: load if available otherwise compute later
flag_l2bounds_compute = 1; % default, will be updated below..
if isempty(filename_l2bounds)
    % assuming Chi squared dist.
    fprintf('\nINFO: l2 bounds will be computed based on the assumption of a chi squared distribution.\n');
else
    fprintf('\nLoading estimates of the l2-bounds ...');
    for iEffCh = 1:nEffectiveChans
        tmpfile=filename_l2bounds(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end));
        if isfile(tmpfile)
            l2EffChansloaded = load(filename_l2bounds(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end)), 'sigmac', 'l2bounds');
            if ~flag_visibility_gridding
                for iCh = 1:numel(effChans2Image{iEffCh})
                    for idSet = 1:nDataSets
                        try % get bounds
                            l2bounds{iEffCh}{iCh}{idSet} = l2EffChansloaded.l2bounds{idSet, iCh};
                            flag_l2bounds_compute = 0;
                            fprintf('\nINFO: l2 bounds loaded successfully.');
                        catch
                            fprintf('\nWARNING: l2 bounds not found. \nChi squared distribution is assumed to determine the l2 bounds.');
                            flag_l2bounds_compute = 1;
                        end
                    end
                end
            else
                try  l2bounds{iEffCh}{1}{1} =  l2EffChansloaded.l2bounds.gridded; % one single block assumed if reduced data
                    fprintf('\nINFO: l2 bounds loaded successfully.\n');
                catch ; fprintf('\nWARNING: l2 bounds not found. \nChi squared distribution is assumed to determine the l2 bounds.\n');
                    flag_l2bounds_compute = 1;
                end
            end
            if ~flag_l2bounds_compute
                try % get noise std
                    if ~flag_visibility_gridding
                        global_sigma_noise(iEffCh, 1) = full(l2EffChansloaded.sigmac);
                    else;  global_sigma_noise(iEffCh, 1) = full(l2EffChansloaded.sigmac.gridded);
                    end
                catch ;  global_sigma_noise(iEffCh, 1) = 1;
                    flag_l2bounds_compute = 1; % assuming a chi squared dist.
                    fprintf('\nWARNING: global sigma noise not found, will assume white data.\n');
                end
                
            end
        else
            fprintf('\nWarning: l2 bounds file not found: %s ',tmpfile);
            flag_l2bounds_compute = 1;
        end
    end
end
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% Auxiliary function needed to select the appropriate workers
% (only needed for 'hs' and 'fhs' algorithms)
switch algo_solver
    case 'sara'
        data_worker_id = @(k) k;
    case 'hs'
        data_worker_id = @(k) k;
    case 'fhs'
        data_worker_id = @(k) k + Qx * Qy;
    otherwise
        error('Undefined algo_solver (available options: ''sara'',''hs'',''fhs'')');
end
%% Get faceting parameter (spectral + spatial)
% fix faceting parameters in case these are not consistent with the
% selected algorithm
if strcmp(algo_solver, 'sara')
    window_type = 'none';
    Qc = 1; % nEffectiveChans;
    Qx = 1;
    Qy = 1;
elseif strcmp(algo_solver, 'hs')
    window_type = 'none';
    Qx = 1;
    Qy = 1;
end
Q = Qx * Qy;

% convert fraction of overlap between consecutive facets into a number of pixels
overlap_size = get_overlap_size([Ny, Nx], [Qy, Qx], overlap_fraction);
if ~strcmp(algo_solver, 'sara')
    fprintf('\nINFO: faceting: number of pixels in overlap: [%s]\n', strjoin(strsplit(num2str(overlap_size))));
end
% index of channels from the subcube to be handled on each data worker
if strcmp(algo_solver, 'sara')
    % for SARA, parallelization is only active over the different wavelet
    % dictionaries composing the SARA prior
    ncores = numel(dict.basis);
    ncores_data = ncores;
else
   freqRangeCores = split_range(ncores_data, nEffectiveChans);
end
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%%  FoV info
param_wproj.FoVx = sin(pixelSize * Nx * pi / 180 / 3600);
param_wproj.FoVy = sin(pixelSize * Ny * pi / 180 / 3600);
param_wproj.uGridSize = 1 / (param_nufft.ox * param_wproj.FoVx);
param_wproj.vGridSize = 1 / (param_nufft.oy * param_wproj.FoVy);
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% setup parpool
try delete(gcp('nocreate'));
end

hpc_cluster = util_set_parpool(algo_solver, ncores_data, Qx * Qy, parcluster);
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% Setup measurement operator and load data
switch algo_solver
    case 'sara'
        for iEffCh = 1:nEffectiveChans

            %% if die/dde calib, get dies/ddes
            if flag_calib.die
                tmpfile = param_preproc.filename_die(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end));
                if isfile(tmpfile)
                    dieloaded = load(param_preproc.filename_die(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end)), 'DIEs');
                else
                    fprintf('\nWarning: DIEs file not found: %s ',tmpfile);
                    dieloaded = [];
                end
            elseif flag_calib.dde
                tmpfile = param_preproc.filename_dde(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end));
                if isfile(tmpfile)
                    ddeloaded = load(param_preproc.filename_dde(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end)), 'DDEs');
                else
                    fprintf('\nWarning: DDEs file not found: %s ',tmpfile);
                    ddeloaded = [];
                end
            end
            %% if data reduction and residual data available, load them
            if flag_visibility_gridding
                try
                    DR_residualdataloaded = load(filename_l2bounds(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end)), 'RESIDUAL_DATA');
                catch ; flag_l2bounds_compute = 1;
                end
            end
            for iCh = 1:numel(effChans2Image{iEffCh})
                for idSet = 1:nDataSets
                    dataloaded = load(dataFilename(idSet, effChans2Image{iEffCh}(iCh)), 'u', 'v', 'w', 'nW', 'y', 'frequency');
                    frequency = dataloaded.frequency;
                    % u v  are in units of the wavelength and will be
                    % normalised between [-pi,pi] for the NUFFT
                    u{iEffCh, 1}{iCh}{idSet} = double(dataloaded.u(:)) * pi / halfSpatialBandwidth;
                    dataloaded.u = [];
                    v{iEffCh, 1}{iCh}{idSet} = -double(dataloaded.v(:)) * pi / halfSpatialBandwidth;
                    dataloaded.v = [];
                    w{iEffCh, 1}{iCh}{idSet} = double(dataloaded.w(:));
                    dataloaded.w = [];
                    nW{iEffCh, 1}{iCh}{idSet} = double(dataloaded.nW(:)); % nW: sqrt(natural weights)
                    dataloaded.nW = [];
                    y{iEffCh, 1}{iCh}{idSet} = double(dataloaded.y(:)) .* nW{iEffCh, 1}{iCh}{idSet}; % data whitening
                    dataloaded.y = [];

                    if flag_calib.die && ~isempty(dieloaded)
                        % if die calib, correct data
                        y{iEffCh, 1}{iCh}{idSet} = y{iEffCh, 1}{iCh}{idSet} ./ dieloaded.DIEs{idSet, iCh};
                        dieloaded.DIEs{idSet, iCh} = [];
                    elseif flag_calib.dde && ~isempty(ddeloaded)
                        % dde calib
                        ddes{iEffCh, 1}{iCh}{idSet} = ddeloaded.DDEs{idSet, iCh};
                        ddeloaded.DDEs{idSet, iCh} = [];
                    end
                    if flag_visibility_gridding
                        % if residual data are available, use them to get
                        % l2 bounds.
                        try preproc_residuals{iEffCh, 1}{iCh}{idSet} = DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh};
                            DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh} = [];
                        catch ;  preproc_residuals = [];
                        end
                    end
                    %% Compute l2bounds if not found
                    if flag_l2bounds_compute
                        global_sigma_noise(iEffCh, 1) = 1;
                        nmeas_curr = numel(y{iEffCh, 1}{iCh}{idSet});
                        l2bounds{iEffCh, 1}{iCh}{idSet} =  sqrt(nmeas_curr + 2 * sqrt(nmeas_curr));
                    end
                end
            end

            %% re-structure data: collapse cells
            u{iEffCh, 1} = vertcat(u{iEffCh, 1}{:});
            v{iEffCh, 1} = vertcat(v{iEffCh, 1}{:});
            w{iEffCh, 1} = vertcat(w{iEffCh, 1}{:});
            nW{iEffCh, 1} = vertcat(nW{iEffCh, 1}{:});
            y{iEffCh, 1} = vertcat(y{iEffCh, 1}{:});
            l2bounds{iEffCh, 1} = vertcat(l2bounds{iEffCh, 1}{:});

            u{iEffCh, 1} = u{iEffCh, 1}(:);
            v{iEffCh, 1} = v{iEffCh, 1}(:);
            w{iEffCh, 1} = w{iEffCh, 1}(:);
            nW{iEffCh, 1} = nW{iEffCh, 1}(:);
            y{iEffCh, 1} =  y{iEffCh, 1}(:);
            l2bounds{iEffCh, 1} = l2bounds{iEffCh, 1}(:);

            if flag_calib.dde
                ddes{iEffCh, 1} = vertcat(ddes{iEffCh, 1}{:});
                ddes{iEffCh, 1} =  ddes{iEffCh, 1}(:);
            else; ddes = [];
            end

            if flag_visibility_gridding
                % if residual data are available
                if exist('preproc_residuals', 'var')
                    preproc_residuals{iEffCh, 1} = vertcat(preproc_residuals{iEffCh, 1}{:});
                    preproc_residuals{iEffCh, 1} = preproc_residuals{iEffCh, 1}(:);
                end
            end
        end

        if flag_visibility_gridding
            % one data block & one channel after visibility gridding
            % ! `G` is the holographic matrix
            % get the reduced data vector `y` and the weighting
            % matrix `Sigma` involved in visibility gridding
            [A, At, G, W, aW, Sigma, y, noise] = util_gen_dr_measurement_operator(y, u, v, w, nW, ...
                1, Nx, Ny, param_nufft, param_wproj, preproc_residuals, ddes);
            % (DR) get the computed l2 bounds from the continuous residual
            % data if available
            if isfield(noise, 'l2bounds') && isfield(noise, 'sigma')
                for iEffCh = 1:nEffectiveChans
                    l2bounds{iEffCh}{1} = noise.l2bounds{iEffCh}{1};
                    global_sigma_noise(iEffCh, 1) = noise.sigma{iEffCh};
                end
            elseif ~isempty(preproc_residuals)
                  fprintf('\nWARNING: computing l2 bounds from the given residual data failed. ');
            end
        else
            [A, At, G, W, aW] = util_gen_measurement_operator(u, v, w, nW, ...
                 1, Nx, Ny, param_nufft, param_wproj, param_precond, ddes);
            Sigma = [];
        end
        clear u v w nW ddes;
        dirty = compute_residual_images(zeros(Ny,Nx), y, A, At, ...
                         G, W, flag_visibility_gridding, Sigma);

    otherwise % 'hs' or 'fhs'
        % create the measurement operator operator in parallel (depending on
        % the algorithm used)
        Sigma = Composite();
        if strcmp(algo_solver, 'hs')
            spmd
                ddes = []; preproc_residuals = [];
                local_fc = (freqRangeCores(labindex, 1):freqRangeCores(labindex, 2));
                effChans2Image_lab = (effChans2Image(local_fc));
                for ifc  = 1:numel(local_fc)
                    %% if die/ddes calib, get dies/ddes                    
                    if flag_calib.die
                        tmpfile =param_preproc.filename_die(effChans2Image{ifc}(1), effChans2Image{ifc}(end));
                        if isfile(tmpfile)
                            dieloaded = load(param_preproc.filename_die(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'DIEs');
                        else
                            fprintf('\nWarning: DIEs file not found: %s ',tmpfile);
                            dieloaded = [];
                        end
                    elseif flag_calib.dde
                        tmpfile =param_preproc.filename_dde(effChans2Image{ifc}(1), effChans2Image{ifc}(end));
                        if isfile(tmpfile)
                            ddeloaded = load(param_preproc.filename_dde(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'DDEs');
                        else
                            fprintf('\nWarning: DDEs file not found: %s ',tmpfile);
                            ddeloaded = [];
                        end
                    end
                    
                    %% if data reduction and residual data available, load them
                    if flag_visibility_gridding
                        try
                            DR_residualdataloaded = load(filename_l2bounds(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'RESIDUAL_DATA');
                        catch ; flag_l2bounds_compute = 1;
                        end
                    end
                    for iCh = 1:numel(effChans2Image_lab{ifc})
                        for idSet = 1:nDataSets
                            %% load data
                            dataloaded = load(dataFilename(idSet, effChans2Image_lab{ifc}(iCh)), 'u', 'v', 'w', 'nW', 'y', 'frequency');
                            % u v  are in units of the wavelength and will be normalised between [-pi,pi] for the NUFFT
                            u{ifc, 1}{iCh}{idSet} = double(dataloaded.u(:)) * pi / halfSpatialBandwidth;
                            dataloaded.u = [];
                            v{ifc, 1}{iCh}{idSet} = -double(dataloaded.v(:)) * pi / halfSpatialBandwidth;
                            dataloaded.v = [];
                            w{ifc, 1}{iCh}{idSet} = double(dataloaded.w(:));
                            dataloaded.w = [];
                            nW{ifc, 1}{iCh}{idSet} = double(dataloaded.nW(:)); % sqrt(natural weights)
                            dataloaded.nW = [];
                            y{ifc, 1}{iCh}{idSet} = double(dataloaded.y(:)) .* nW{ifc, 1}{iCh}{idSet}; % data whitening
                            dataloaded.y = [];
                            % Output of the pre-processing step if any
                            if flag_calib.die && ~isempty(dieloaded)
                                % if die calib, correct data
                                y{ifc, 1}{iCh}{idSet} = y{ifc, 1}{iCh}{idSet} ./ dieloaded.DIEs{idSet, iCh};
                                dieloaded.DIEs{idSet, iCh} = [];
                            elseif flag_calib.dde  && ~isempty(ddeloaded)
                                % if dde calib, get ddes
                                ddes{ifc, 1}{iCh}{idSet} = ddeloaded.DDEs{idSet, iCh};
                                ddeloaded.DDEs{idSet, iCh} = [];
                            end
                            if flag_visibility_gridding
                                % if residual data are available, use them to get
                                % l2 bounds.
                                try preproc_residuals{ifc, 1}{iCh}{idSet} = DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh};
                                    DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh} = [];
                                catch ; preproc_residuals = [];
                                end
                            end

                            %% Compute l2bounds if not found
                            if flag_l2bounds_compute
                                global_sigma_noise_cmpst = 1;
                                nmeas_curr = numel(y{ifc, 1}{iCh}{idSet});
                                l2bounds{ifc, 1}{iCh}{idSet} =  sqrt(nmeas_curr + 2 * sqrt(nmeas_curr));
                            end
                        end
                    end
                    %% re-structure data: collapse cells
                    u{ifc, 1} = vertcat(u{ifc, 1}{:});
                    v{ifc, 1} = vertcat(v{ifc, 1}{:});
                    w{ifc, 1} = vertcat(w{ifc, 1}{:});
                    nW{ifc, 1} = vertcat(nW{ifc, 1}{:});
                    y{ifc, 1} = vertcat(y{ifc, 1}{:});
                    l2bounds{ifc, 1} = vertcat(l2bounds{ifc, 1}{:});
                    %
                    u{ifc, 1} = u{ifc, 1}(:);
                    v{ifc, 1} = v{ifc, 1}(:);
                    w{ifc, 1} = w{ifc, 1}(:);
                    nW{ifc, 1} = nW{ifc, 1}(:);
                    y{ifc, 1} =  y{ifc, 1}(:);
                    l2bounds{ifc, 1} = l2bounds{ifc, 1}(:);

                    if flag_calib.dde
                        ddes{ifc, 1} = vertcat(ddes{ifc, 1}{:});
                        ddes{ifc, 1} = ddes{ifc, 1}(:);
                    else; ddes = [];
                    end

                    if flag_visibility_gridding
                        try preproc_residuals{ifc, 1} = vertcat(preproc_residuals{ifc, 1}{:});
                            preproc_residuals{ifc, 1} = preproc_residuals{ifc, 1}(:);
                        catch ; preproc_residuals = [];
                        end
                    end
                end

                if flag_visibility_gridding
                    % one data block & one channel after visibility gridding
                    % ! `G` is the holographic matrix
                    % get the reduced data vector `y` and the weighting
                    % matrix `Sigma` involved in visibility gridding
                    fprintf('\nCompute the holographic matrix  .. ');
                    [A, At, G, W, aW, Sigma, y] = util_gen_dr_measurement_operator(y, u, v, w, nW, ...
                         numel(local_fc), Nx, Ny, param_nufft, param_wproj, preproc_residuals, ddes);
                else
                    % ! `G` is the degridding matrix
                    [A, At, G, W, aW] = util_gen_measurement_operator(u, v, w, nW, ...
                         numel(local_fc), Nx, Ny, param_nufft, param_wproj, param_precond, ddes);
                    Sigma = [];
                end;   u = []; v = []; w = []; nW = []; ddes = [];
                dirty = compute_residual_images(squeeze(zeros(Ny,Nx,numel(y))), y, A, At, ...
                                                G, W, flag_visibility_gridding, Sigma);
            end
        else
            Sigma = Composite();
            spmd
                ddes = []; preproc_residuals = [];
                % define operator on data workers only
                if labindex > Q
                    local_fc = (freqRangeCores(labindex - Q, 1):freqRangeCores(labindex - Q, 2));
                    effChans2Image_lab = (effChans2Image(local_fc));
                    for ifc  = 1:numel(local_fc)
                        %% if die/ddes calib, get dies/ddes
                        if flag_calib.die
                            tmpfile =param_preproc.filename_die(effChans2Image{ifc}(1), effChans2Image{ifc}(end));
                            if isfile(tmpfile)
                                dieloaded = load(param_preproc.filename_die(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'DIEs');
                            else
                                fprintf('\nWarning: DIEs file not found: %s ',tmpfile);
                                dieloaded = [];
                            end
                        elseif flag_calib.dde
                            tmpfile =param_preproc.filename_dde(effChans2Image{ifc}(1), effChans2Image{ifc}(end));
                            if isfile(tmpfile)
                                ddeloaded = load(param_preproc.filename_dde(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'DDEs');
                            else
                                fprintf('\nWarning: DDEs file not found: %s ',tmpfile);
                                ddeloaded = [];
                            end
                        end                        
                        %% if data reduction and residual data available, load them
                        if flag_visibility_gridding
                            try
                                DR_residualdataloaded = load(filename_l2bounds(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'RESIDUAL_DATA');
                            catch ; flag_l2bounds_compute = 1;
                            end
                        end
                        for iCh = 1:numel(effChans2Image_lab{ifc})
                            for idSet = 1:nDataSets
                                %% load data
                                dataloaded = load(dataFilename(idSet, effChans2Image_lab{ifc}(iCh)), 'u', 'v', 'w', 'nW', 'y', 'frequency');
                                % u v  are in units of the wavelength and will be normalised between [-pi,pi] for the NUFFT
                                u{ifc, 1}{iCh}{idSet} = double(dataloaded.u(:)) * pi / halfSpatialBandwidth;
                                dataloaded.u = [];
                                v{ifc, 1}{iCh}{idSet} = -double(dataloaded.v(:)) * pi / halfSpatialBandwidth;
                                dataloaded.v = [];
                                w{ifc, 1}{iCh}{idSet} = double(dataloaded.w(:));
                                dataloaded.w = [];
                                nW{ifc, 1}{iCh}{idSet} = double(dataloaded.nW(:)); % sqrt(natural weights)
                                dataloaded.nW = [];
                                y{ifc, 1}{iCh}{idSet} = double(dataloaded.y(:)) .* nW{ifc, 1}{iCh}{idSet}; % data whitening
                                dataloaded.y = [];
                                % Output of the pre-processing step if any
                                if flag_calib.die && ~isempty(dieloaded)
                                    % if die calib, correct data
                                    y{ifc, 1}{iCh}{idSet} = y{ifc, 1}{iCh}{idSet} ./ dieloaded.DIEs{idSet, iCh};
                                    dieloaded.DIEs{idSet, iCh} = [];
                                elseif flag_calib.dde && ~isempty(ddeloaded)
                                    % if dde calib, get ddes
                                    ddes{ifc, 1}{iCh}{idSet} = ddeloaded.DDEs{idSet, iCh};
                                    ddeloaded.DDEs{idSet, iCh} = [];
                                end
                                if flag_visibility_gridding
                                    % if residual data are available, use them to get
                                    % l2 bounds.
                                    try preproc_residuals{ifc, 1}{iCh}{idSet} = DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh};
                                        DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh} = [];
                                    catch ; preproc_residuals = [];
                                    end
                                end

                                %% Compute l2bounds if not found
                                if flag_l2bounds_compute
                                    global_sigma_noise_cmpst = 1;
                                    nmeas_curr = numel(y{ifc, 1}{iCh}{idSet});
                                    l2bounds{ifc, 1}{iCh}{idSet} =  sqrt(nmeas_curr + 2 * sqrt(nmeas_curr));
                                end
                            end
                        end
                        %% re-structure data: collapse cells
                        u{ifc, 1} = vertcat(u{ifc, 1}{:});
                        v{ifc, 1} = vertcat(v{ifc, 1}{:});
                        w{ifc, 1} = vertcat(w{ifc, 1}{:});
                        nW{ifc, 1} = vertcat(nW{ifc, 1}{:});
                        y{ifc, 1} = vertcat(y{ifc, 1}{:});
                        l2bounds{ifc, 1} = vertcat(l2bounds{ifc, 1}{:});
                        %
                        u{ifc, 1} = u{ifc, 1}(:);
                        v{ifc, 1} = v{ifc, 1}(:);
                        w{ifc, 1} = w{ifc, 1}(:);
                        nW{ifc, 1} = nW{ifc, 1}(:);
                        y{ifc, 1} =  y{ifc, 1}(:);
                        l2bounds{ifc, 1} = l2bounds{ifc, 1}(:);

                        if flag_calib.dde 
                            ddes{ifc, 1} = vertcat(ddes{ifc, 1}{:});
                            ddes{ifc, 1} = ddes{ifc, 1}(:);
                        else; ddes = [];
                        end

                        if flag_visibility_gridding
                            try preproc_residuals{ifc, 1} = vertcat(preproc_residuals{ifc, 1}{:});
                                preproc_residuals{ifc, 1} = preproc_residuals{ifc, 1}(:);
                            catch ; preproc_residuals = [];
                            end
                        end
                    end

                    if flag_visibility_gridding
                        % one data block & one channel after visibility gridding
                        % ! `G` is the holographic matrix
                        % get the reduced data vector `y` and the weighting
                        % matrix `Sigma` involved in visibility gridding
                        fprintf('\nCompute the holographic matrix  ..');
                        [A, At, G, W, aW, Sigma, y] = util_gen_dr_measurement_operator(y, u, v, w, nW, ...
                            numel(local_fc), Nx, Ny, param_nufft, param_wproj,  preproc_residuals, ddes);
                    else
                        % ! `G` is the degridding matrix
                        [A, At, G, W, aW] = util_gen_measurement_operator(u, v, w, nW, ...
                            numel(local_fc), Nx, Ny, param_nufft, param_wproj,  param_precond, ddes);
                        Sigma = [];
                    end; u = []; v = []; w = []; nW = []; ddes = [];
                    
                    dirty = compute_residual_images(squeeze(zeros(Ny,Nx,numel(y))), y, A, At, ...
                                                G, W, flag_visibility_gridding, Sigma);
                end
            end
        end; clear local_fc  u v w nW dataSpWinloaded ddes;
end
% Free memory
clear  param_wproj param_preproc  preproc_residuals;
  
fprintf('\nData loaded successfully.');
try
    subcube_channels = 1:nEffectiveChans;
    if strcmp(algo_solver, 'sara')
        dirty_file= fullfile(results_path, ...
            strcat('dirty_', algo_solver, ...
            '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
            '_chs', num2str(effChans2Image{1}(1)), '-', num2str(effChans2Image{1}(end)), '.fits'));
        fitswrite(dirty,dirty_file)
    else
        for k = 1:ncores_data
            dirty_file = fullfile(results_path, strcat('dirty_', ...
                '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
                '_ch', num2str(k),'.fits'));
            dirty_tmp = dirty{data_worker_id(k)};
            fitswrite(dirty_tmp,dirty_file)
        end
    end
    fprintf('\nDirty images saved.')
catch, fprintf('\nWARNING: error encountered when saving dirty images.');
end
%% load l2 bounds (generate only full spectral dataset)

switch algo_solver
    case 'sara'
        % ! to be verified
        % all the variables are stored on the main process for sara
        epsilons = l2bounds(1);
    otherwise
        epsilons = Composite();
        for k = 1:ncores_data
            epsilons{data_worker_id(k)} = l2bounds{data_worker_id(k)};
            if exist('global_sigma_noise_cmpst', 'var')
                if k == 1; global_sigma_noise = global_sigma_noise_cmpst{data_worker_id(k)};
                else; global_sigma_noise = [global_sigma_noise; global_sigma_noise_cmpst{data_worker_id(k)}];
                end
            end
        end; clear global_sigma_noise_cmpst;
end

%% Compute operator norm
subcube_channels = 1:nEffectiveChans;
if strcmp(algo_solver, 'sara')
    if flag_computeOperatorNorm
        fprintf('\nINFO: computing operator''s spectral norm .. ')
        tic
        [Anorm, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-6, 200, flag_visibility_gridding, Sigma); % AD: changed tolerance to 1e-6 instead of 1e-8
        save(fullfile(results_path, ...
            strcat('Anorm_', algo_solver, ...
            '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
            '_chs', num2str(effChans2Image{1}(1)), '-', num2str(effChans2Image{1}(end)), '.mat')), ...
            '-v7.3', 'Anorm', 'squared_operator_norm', 'rel_var', ...
            'squared_operator_norm_precond', 'rel_var_precond');
        clear rel_var;
        fprintf('done in %d sec.\n',floor(toc))
    else
        load(fullfile(results_path, ...
            strcat('Anorm_', algo_solver, ...
            '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
            '_chs', num2str(effChans2Image{1}(1)), '-', num2str(effChans2Image{1}(end)), '.mat')), ...
            'Anorm', 'squared_operator_norm_precond', 'squared_operator_norm');
    end
else
    if flag_computeOperatorNorm
        fprintf('\nINFO: computing operator''s spectral norm .. ')
        tic
        spmd
            if labindex > Qx * Qy * strcmp(algo_solver, 'fhs')
                [An, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-8, 200, flag_visibility_gridding, Sigma);
            end
        end
        % save operator norm from the different subcubes into a single .mat
        % file
        opnormfile = matfile(fullfile(results_path, strcat('Anorm', ...
            '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
            '_L', num2str(nEffectiveChans), '.mat')), 'Writable', true);

        opnormfile.squared_operator_norm = zeros(nEffectiveChans, 1);
        opnormfile.rel_var = zeros(nEffectiveChans, 1);
        opnormfile.squared_operator_norm_precond = zeros(nEffectiveChans, 1);
        opnormfile.rel_var_precond = zeros(nEffectiveChans, 1);

        Anorm = 0;
        for k = 1:ncores_data
            opnormfile.squared_operator_norm(freqRangeCores(k, 1):freqRangeCores(k, 2), 1) = squared_operator_norm{data_worker_id(k)};
            opnormfile.rel_var(freqRangeCores(k, 1):freqRangeCores(k, 2), 1) = rel_var{data_worker_id(k)};

            opnormfile.squared_operator_norm_precond(freqRangeCores(k, 1):freqRangeCores(k, 2), 1) = squared_operator_norm_precond{data_worker_id(k)};
            opnormfile.rel_var_precond(freqRangeCores(k, 1):freqRangeCores(k, 2), 1) = rel_var_precond{data_worker_id(k)};

            Anorm = max(Anorm, An{data_worker_id(k)});
        end
        clear An rel_var rel_var_precond squared_operator_norm_precond;
        squared_operator_norm =   opnormfile.squared_operator_norm(subcube_channels, 1);
        squared_operator_norm_precond = opnormfile.squared_operator_norm_precond(subcube_channels, 1);

        fprintf('done in %d sec.\n',floor(toc))
    else
        opnormfile = matfile(fullfile(results_path, strcat('Anorm', ...
            '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
            '_L', num2str(nEffectiveChans), '.mat')));

        squared_operator_norm_precond = opnormfile.squared_operator_norm_precond(subcube_channels, 1);
        rel_var_precond = opnormfile.rel_var_precond(subcube_channels, 1);
        Anorm = max(squared_operator_norm_precond .* (1 + rel_var_precond));
        squared_operator_norm = opnormfile.squared_operator_norm(subcube_channels, 1);
    end
end

fprintf('INFO: Convergence parameter (measurement operator''s norm squared): %e \n', Anorm);
%% Setup name of results file

if strcmp(algo_solver, 'sara') % AD
    temp_results_name = @(nEffectiveChans) strcat(srcName, '_', ...
        algo_solver, '_', num2str(pixelSize), 'asec', ...
        '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
        '_ind', num2str(subcubeInd), '_chs', num2str(effChans2Image{1}(1)), '-', num2str(effChans2Image{1}(end)), ...
        '_g', num2str(gam), '_gb', num2str(gam_bar), '_');
else
    temp_results_name = @(nEffectiveChans) strcat(srcName, '_', ...
        algo_solver, '_', num2str(pixelSize), 'asec', ...
        '_Ny', num2str(Ny), '_Nx', num2str(Nx), '_L', num2str(nEffectiveChans), ...
        '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
        '_ind', num2str(subcubeInd), '_g', num2str(gam), '_gb', num2str(gam_bar), ...
        '_overlap', strjoin(strsplit(num2str(overlap_fraction)), '_'));
end

warm_start = @(nEffectiveChans) strcat(temp_results_name(nEffectiveChans), '_rw', num2str(flag_reweighting), '.mat');
%% Regularization parameters and solver

% estimate noise level (set regularization parameters to the same value)
% compute sig and sig_bar (estimate of the "noise level" in "SVD" and
% SARA space) involved in the reweighting scheme

if strcmp(algo_solver, 'sara')
    % SARA dicionary (created out of the solver for SARA)
    dwtmode('zpd');
    [Psi1, Psit1] = op_p_sp_wlt_basis(dict.basis, dict.nlevel, Ny, Nx);
    P = numel(Psi1);
    Psi = cell(P, 1);
    Psit = cell(P, 1);
    s = zeros(P, 1); % number of wavelet coefficients for each dictionary
    for k = 1:P
        f = '@(x_wave) HS_forward_sparsity(x_wave,Psi1{';
        f = sprintf('%s%i},Ny,Nx);', f, k);
        Psi{k} = eval(f);
        s(k) = size(Psit1{k}(zeros(Ny, Nx, 1)), 1);
        ft = ['@(x) HS_adjoint_sparsity(x,Psit1{' num2str(k) '},s(' num2str(k) '));'];
        Psit{k} = eval(ft);
    end
    % noise level / regularization parameter
    sig = compute_noise_level_sara(global_sigma_noise, squared_operator_norm);
    % apply multiplicative factor for the regularization parameter (if needed)
    mu = gam * sig;
    fprintf('Noise level: sig = %e\n', sig);
    fprintf('Additional multiplicative regularisation factor gam = %e\n', gam);
    fprintf('Regularization parameter mu = %e\n', mu);
    fprintf('Algo: %s, alpha = %.4e, mu = %.4e, sig = %.4e\n', algo_solver, gam, mu, sig);
end

if strcmp(algo_solver, 'hs') || strcmp(algo_solver, 'fhs')
    % noise level / regularization parameter
    [sig, sig_bar, mu_chi, sig_chi, sig_sara] = ...
        compute_noise_level(Ny, Nx, nEffectiveChans, global_sigma_noise(:), ...
        algo_solver, Qx, Qy, overlap_size, squared_operator_norm(:));
    % apply multiplicative factor for the regularization parameters (if needed)
    mu_bar = gam_bar * sig_bar;
    mu = gam * sig;
    fprintf('mu_chi = %.4e, sig_chi = %.4e, sig_sara = %.4e\n', mu_chi, sig_chi, sig_sara);
    fprintf('Noise levels: sig = %.4e, sig_bar = [%.4e, %.4e]\n', sig, min(sig_bar), max(sig_bar));
    fprintf('Additional multiplicative actors gam = %.4e, gam_bar = %.4e\n', gam, gam_bar);
    fprintf('Regularization parameters: mu = %.4e, mu_bar = %.4e\n', mu, mu_bar);
    fprintf('Algo: %s, gam = %.4e, gam_bar = %.4e, mu = %.4e, mu_bar = [%.4e, %.4e]\n', algo_solver, gam, gam_bar, mu, min(mu_bar), max(mu_bar));
end

% moved here for transprancy: define parameters for the solver
% List of default solver-specific parameters (reweighting, pdfb, ellispsoid
% prjection, epsilon update scheme).
% * general
% estimate of the noise level in SARA space
param_solver.reweighting_sig = sig;
if ~strcmp(algo_solver, 'sara')
    param_solver.reweighting_sig_bar = sig_bar; % estimate of the noise level in "SVD" spaces
end
param_solver.nu0 = 1; % bound on the norm of the Identity operator
param_solver.nu1 = 1; % bound on the norm of the operator Psi
param_solver.nu2 = squared_operator_norm_precond; % upper bound on the norm of the measurement operator
if ~strcmp(algo_solver, 'sara')
    param_solver.gamma0 = mu_bar; % regularization parameter nuclear norm
end
param_solver.gamma = mu; % regularization parameter l21-norm (soft th parameter) ! for SARA, take the value given as an input to the solver
param_solver.cube_id = subcubeInd; % id of the cube to be reconstructed
param_solver.backup_frequency = 1; % PA: checkpoint frequency! AD :????

param_solver.alph = gam;
param_solver.alph_bar = gam_bar;
% temp filenames
name_checkpoint = fullfile(auxiliary_path, temp_results_name(nEffectiveChans));
name_warmstart = fullfile(auxiliary_path, warm_start(nEffectiveChans));
%% Solver
if flag_solveMinimization
    if strcmp(algo_solver, 'sara')
        disp('SARA');
        disp('-----------------------------------------');
        % ! in this case, ncores_data corresponds
        % ! to the number of workers for the wavelet transform (9 maximum)
        xsol = sara(y, epsilons, A, At, aW, G, W, Psi, Psit, ...
            param_solver, name_warmstart, name_checkpoint, gam, ...
            flag_visibility_gridding, Sigma, xinit);
        fitswrite(xsol, fullfile(auxiliary_path, strcat('x_', srcName, '_', algo_solver, ...
            '_', num2str(pixelSize), 'asec', ...  % '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
            '_chs', num2str(effChans2Image{1}(1)), '-', num2str(effChans2Image{1}(end)), ...
            '_gam', num2str(gam), ...
            '.fits')));
    else
        %%
        % spectral tesselation (non-overlapping)
        channel_chunks = cell(ncores_data, 1);
        for k = 1:ncores_data
            channel_chunks{k} = freqRangeCores(k, 1):freqRangeCores(k, 2);
        end

        %
        xinit = reshape(xinit, ImageCubeDims);
        switch algo_solver
            case 'hs'
                disp('HyperSARA');
                disp('-----------------------------------------');
                xsol = hyperSARA(y, epsilons, ...
                    A, At, aW, G, W, param_solver, ...
                    ncores_data, dict.basis, dict.nlevel, channel_chunks, ...
                    nEffectiveChans, Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flag_visibility_gridding, Sigma, ...
                    xinit);
            case 'fhs'
                disp('Faceted HyperSARA');
                disp('-----------------------------------------');

                xsol = facetHyperSARA(y, epsilons, ...
                    A, At, aW, G, W, param_solver, Qx, Qy, ncores_data, ...
                    dict.basis, dict.filter_length, dict.nlevel, window_type, ...
                    channel_chunks, nEffectiveChans, overlap_size, gam, gam_bar, ...
                    Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flag_visibility_gridding, Sigma, ...
                    xinit);
            otherwise; error('Unknown solver version.');
        end
        fitswrite(xsol, fullfile(auxiliary_path, strcat('x_', srcName, '_', algo_solver, ...
            '_', num2str(pixelSize), 'asec', ...
            '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
            '_ind', num2str(subcubeInd), ...
            '_gam', num2str(gam), '_gambar', num2str(gam_bar), ...
            '_overlap', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
            '.fits')));
    end
end
