function imaging(srcName, datasetsNames, dataFilename, ...
    effChans2Image, param_solver, param_global, ...
    param_nufft, param_precond, param_wproj, dict)
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
%     Names of the different datasets to be imaged (e.g., in the cases of
%     different acquisition configurations or different observation times),
%     to be set to ``[]`` if one data set is imaged and no name tag of
%     the associated ``MS`` is given during data extraction.
% dataFilenames: cell of string
%     Function handle taking a channel index as an input and returning the
%     name of the associated data ``.mat`` file (one file per frequency channel
%     per dataset), following the nomenclature adopted during data extraction.
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
% param_global.adjust_flag_noise : bool
%     Flag to activate adaptive estimation of the noise level & paramerter
%     adjustment.
% param_global.measop_flag_visibility_gridding : bool
%     Flag to activate data dimensionality reduction via visibility gridding in the
%     definition of the measurement operator.
% param_global.data_flag_apply_imaging_weights : bool
%     Flag to read and apply imaging weights (e.g. Briggs, uniform) to data.
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
% param_global.reg_nReweights: int
%     Maximum number of reweight.
% param_global.algo_flag_computeOperatorNorm : bool
%     Flag to trigger the computation of the norm of the (preconditionned) measurement
%     operator. Default set to ``true``. If set to ``false``,
%     MATLAB will look for a file where this
%     quantity has been saved (save and computation steps are triggered in
%     :mat:func:`imaging.imaging`).
% param_global.algo_flag_saveOperatorNorm : bool
%     Flag to save the norm of the (preconditionned) measurement
%     operator. Default set to ``false``. If set to ``true``,
%     MATLAB will save the quantity in a ``.mat`` file.
% param_global.algo_flag_solveMinimization : bool
%     Flag triggering the solver (``"fhs"``, ``"hs"`` or ``"sara"``).
% param_global.ncores_data : int
%     Number of cores handlig the data fidelity terms (data cores). For
%     Faceted HyperSARA, the total number of cores used
%     is ``Qx*Qy + ncores_data + 1``. For SARA and HyperSARA, represents
%     the number of cores used for the parallelization.
% param_global.parcluster : string
%     Name of the parallel parcluster profile to launch the parpool.
%     By default ``"local"`` profile  is used. The user should set it to
%     the name of the slurm parcluster profile created on his/her HPC
%     machine, if prefered.
% param_global.preproc_filename_noise_std : anonymous function
%     Function handle, taking the indices of the physical channel, and the
%     and the dataset, returning a string corresponding to the name of a file
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
% Note
% ----
% - All the results will be stored in a directory of th
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
default_superresolution = 2;
%% Paths
format compact;
% Project dir.
if ~isfield(param_global, 'main_dir'); param_global.main_dir = [pwd, filesep]; end
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
addpath([project_dir, filesep, 'lib', filesep, 'RI-measurement-operator', filesep, 'nufft']);
addpath([project_dir, filesep, 'lib', filesep, 'RI-measurement-operator', filesep, 'lib', filesep, 'utils']);
addpath([project_dir, filesep, 'lib', filesep, 'RI-measurement-operator', filesep, 'lib', filesep, 'operators']);
addpath([project_dir, filesep, 'lib', filesep, 'RI-measurement-operator', filesep, 'lib', filesep, 'ddes_utils']);
addpath([project_dir, filesep, 'lib', filesep, 'SARA-dictionary', filesep, 'src']);
addpath([project_dir, filesep, 'src']);
addpath([project_dir, filesep, 'src', filesep, 'heuristics', filesep]);
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% check input params
% Image resolution & dimensions
% ! parameters and default values (keep only what is required)
if ~isfield(param_global, 'im_pixelSize');  param_global.im_pixelSize = []; end

% Prior: sparsity dict.
if ~isfield(dict, 'nlevel'); dict.nlevel = 4;  end
if ~isfield(dict, 'basis')
    dict.basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};
    dict.filter_length = [2, 4, 6, 8, 10, 12, 14, 16, 0];
end

% Prior: wideband facet related
max_filter_length = max(dict.filter_length(:));
if ~isfield(param_global, 'facet_subcubeInd');  param_global.facet_subcubeInd = 0; end
if ~isfield(param_global, 'facet_Qx')
    param_global.facet_Qx = sdwt2_max_nfacets(param_global.im_Nx, ndict.level, max_filter_length);
else
    nfacets_max = sdwt2_test_nfacets(param_global.im_Nx, param_global.facet_Qx, dict.nlevel, max_filter_length);
    if nfacets_max < param_global.facet_Qx
        param_global.facet_Qx = nfacets_max;
        fprintf('\nWARNING: facet dimension is too small, max. admissible value used: param_global.facet_Qx=%d ', param_global.facet_Qx);
    end
end
if ~isfield(param_global, 'facet_Qy')
    param_global.facet_Qy = sdwt2_max_nfacets(param_global.im_Ny, ndict.level, max_filter_length);
else
    nfacets_max = sdwt2_test_nfacets(param_global.im_Ny, param_global.facet_Qy, dict.nlevel, max_filter_length);
    if nfacets_max < param_global.facet_Qy
        param_global.facet_Qy = nfacets_max;
        fprintf('\nWARNING: facet dimension is too small, max. admissible value used: param_global.facet_Qy=%d ', param_global.facet_Qy);
    end
end
if ~isfield(param_global, 'facet_Qc'); param_global.facet_Qc = 1; end
if ~isfield(param_global, 'facet_window_type'); param_global.facet_window_type = 'triangular'; end
if ~isfield(param_global, 'facet_overlap_fraction'); param_global.facet_overlap_fraction = [0.2, 0.2]; end

% Prior: reg params
if ~isfield(param_global, 'reg_gam'); param_global.reg_gam = 1;  end
if ~isfield(param_global, 'reg_gam_bar'); param_global.reg_gam_bar = 1;  end
if ~isfield(param_global, 'reg_flag_reweighting'); param_global.reg_flag_reweighting = true; end

if param_global.reg_flag_reweighting &&  ~isfield(param_global, 'reg_nReweights') && ~isfield(param_solver, 'reweighting_max_iter')
    param_solver.reweighting_max_iter = 10;
elseif isfield(param_global, 'reg_nReweights')
    % overwrite the default value set in .json file
    param_solver.reweighting_max_iter = param_global.reg_nReweights;
end

% Noise level adjustement
if ~isfield(param_global, 'adjust_flag_noise'); param_global.adjust_flag_noise = false; end

% Algo
if numel(effChans2Image) == 1
    param_global.algo_solver = 'sara';
    fprintf('\nMonochromatic imaging: `SARA` approach is used for imaging.');
elseif ~isfield(param_global, 'algo_solver')
    param_global.algo_solver = 'fhs';
    fprintf('\nWideband imaging: `Faceted HyperSARA` approach is used for imaging.');
elseif numel(effChans2Image) > 1 && strcmp(param_global.algo_solver, 'sara')
    param_global.algo_solver = 'fhs';
    fprintf('\nWideband imaging: `Faceted HyperSARA` approach is used for imaging.');
end
if ~isfield(param_global, 'algo_flag_solveMinimization'); param_global.algo_flag_solveMinimization = true; end
if ~isfield(param_global, 'algo_flag_computeOperatorNorm'); param_global.algo_flag_computeOperatorNorm = true; end
if ~isfield(param_global, 'algo_flag_saveOperatorNorm'); param_global.algo_flag_saveOperatorNorm = false; end

if ~isfield(param_global, 'algo_ncores_data'); param_global.algo_ncores_data = numel(effChans2Image);  end

% Measurement operator: visibility gridding
if ~isfield(param_global, 'measop_flag_visibility_gridding'); param_global.measop_flag_visibility_gridding = false; end

% data and measurement operator : imaging weights
if ~isfield(param_global, 'data_flag_apply_imaging_weights'); param_global.data_flag_apply_imaging_weights = false; end

% Measurement operator: W-projection
if ~isfield(param_wproj, 'measop_flag_wproj'); param_wproj.measop_flag_wproj = false; end
if param_wproj.measop_flag_wproj
    if ~isfield(param_wproj, 'CEnergyL2'); param_wproj.CEnergyL2 = 1 - 1e-4; end
    if ~isfield(param_wproj, 'GEnergyL2'); param_wproj.GEnergyL2 = 1 - 1e-4; end
else
    param_wproj.CEnergyL2 = 1;
    param_wproj.GEnergyL2 = 1;
end

% Input filenames from pre-processing
if ~isfield(param_global, 'preproc_filename_model'); param_global.preproc_filename_model = []; end
if ~isfield(param_global, 'preproc_filename_noise_std'); param_global.preproc_filename_noise_std = []; end
if ~isfield(param_global, 'preproc_filename_cal_solutions'); param_global.preproc_filename_cal_solutions = []; end

% Output filenames
if ~isfield(param_global, 'exp_type'); param_global.exp_type = '';  end %

% name of the parallel parcluster profile
if ~isfield(param_global, 'parcluster'); param_global.parcluster = 'local';
end
%% rest of the paths -----------------------------------------------------------------%
% ------------------------------------------------------------------------%
algo_solver = param_global.algo_solver;
switch algo_solver
    case 'sara'
        addpath([project_dir, filesep, 'src', filesep, 'sara']);
    case 'hs'
        addpath([project_dir, filesep, 'src', filesep, 'hs', filesep]);
    case 'fhs'
        addpath([project_dir, filesep, 'src', filesep, 'fhs', filesep]);
end
% setting paths to results and reference image cube
results_path = fullfile(project_dir, 'results', srcName);
mkdir(results_path);
auxiliary_path = fullfile(results_path, algo_solver);
mkdir(auxiliary_path);
fprintf('\nINFO: results will be saved in: \n %s .\n', auxiliary_path);

% -------------------------------------------------------------------------%
%% get params
% image dimensions & resolution
Nx = param_global.im_Nx;
Ny = param_global.im_Ny;
pixelSize  = param_global.im_pixelSize;
if isempty(pixelSize); imResolution = 'nominal';
else; imResolution = 'user_defined';
end
switch imResolution
    case 'nominal'
        fprintf('\nWARNING: pixelsize not found, default is considered.');
    otherwise
        fprintf('\nINFO: pixelsize: %f arcsec.\n', pixelSize);
end
% Prior: wideband facet related
subcubeInd = param_global.facet_subcubeInd;
% Prior: faceting: de-activate for SARA & HyperSARA
% fix faceting parameters in case these are not consistent with the
% selected algorithm
switch algo_solver

    case 'fhs'
        Qx = param_global.facet_Qx;
        Qy = param_global.facet_Qy;
        Qc = param_global.facet_Qc;
        window_type = param_global.facet_window_type;
        overlap_fraction = param_global.facet_overlap_fraction;
    otherwise
        window_type = 'none';
        Qc = 1;
        Qx = 1;
        Qy = 1;
        overlap_fraction = [0 0];
end

% Prior: reg params
gam = param_global.reg_gam;
gam_bar =  param_global.reg_gam_bar;
flag_reweighting = param_global.reg_flag_reweighting;
if flag_reweighting;  reg_nReweights = param_global.reg_nReweights;
end

% Noise estimation
adjust_flag_noise = param_global.adjust_flag_noise;

% Measurement operator settings
flag_visibility_gridding = param_global.measop_flag_visibility_gridding; % data reduction
% data and measurement operator
flag_apply_imaging_weights = param_global.data_flag_apply_imaging_weights;
% Algo
ncores_data = param_global.algo_ncores_data;
flag_computeOperatorNorm = param_global.algo_flag_computeOperatorNorm;
flag_saveOperatorNorm = param_global.algo_flag_saveOperatorNorm;
flag_solveMinimization =  param_global.algo_flag_solveMinimization;
% Pre-processing step
param_preproc_files.dde  = param_global.preproc_filename_cal_solutions;
param_preproc_files.noise_std = param_global.preproc_filename_noise_std;
filename_init_model = param_global.preproc_filename_model;

% parcluster
parcluster = param_global.parcluster;
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% info
disp('Imaging configuration');
disp(['Source name: ', srcName]);
disp(['Image cube size: ', num2str(Ny), ' x ', num2str(Nx), ' x ', num2str(numel(effChans2Image))]);
disp(['Algorithm version: ', algo_solver]);

if ~strcmp(algo_solver, 'sara')
    disp(['Number of facets Qy x Qx : ', num2str(Qy), ' x ', num2str(Qx)]);
    disp(['Number of spectral facets Qc : ', num2str(Qc)]);
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
    fprintf('\nINFO: default pixelsize: %f arcsec, that is %d x nominal resolution at the highest freq.', pixelSize, default_superresolution);
end
halfSpatialBandwidth = (180 / pi) * 3600 / (pixelSize) / 2;
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% check noise estimation
if isempty(param_preproc_files.noise_std)
    flag_noise_theoretical = 1;
    fprintf('\nINFO: theoretical noise assumed.\n');
else; flag_noise_theoretical = 0;
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
%% Get faceting parameter (spectral + spatial) & index of channels from the subcube to be handled on each data worker

switch algo_solver
    case 'sara'
        % for SARA, parallelization is only active over the different wavelet
        % dictionaries composing the SARA prior
        ncores = numel(dict.basis);
        ncores_data = ncores;
    otherwise
        freqRangeCores = split_range(ncores_data, nEffectiveChans);
        % convert fraction of overlap between consecutive facets into a number of pixels
        overlap_size = get_overlap_size([Ny, Nx], [Qy, Qx], overlap_fraction);
        fprintf('\nINFO: faceting: number of pixels in overlap: [%s]\n', strjoin(strsplit(num2str(overlap_size))));

end
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%%  FoV info for w-proj
param_wproj.FoVx = sin(pixelSize * Nx * pi / 180 / 3600);
param_wproj.FoVy = sin(pixelSize * Ny * pi / 180 / 3600);
param_wproj.uGridSize = 1 / (param_nufft.ox * param_wproj.FoVx);
param_wproj.vGridSize = 1 / (param_nufft.oy * param_wproj.FoVy);
param_wproj.halfSpatialBandwidth = halfSpatialBandwidth;
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% Preconditionning params
param_precond.N = param_global.im_Nx * param_global.im_Ny; % number of Fourier points (oversampled plane)
param_precond.Nox = param_nufft.ox * param_global.im_Nx;
param_precond.Noy = param_nufft.oy * param_global.im_Ny; %% Setup measurement operator and load data

%% Measurement operator
if isempty(flag_apply_imaging_weights); flag_apply_imaging_weights = false;
end
param_nufft.flag_apply_imaging_weights = flag_apply_imaging_weights;

% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% Building the measurement operator
switch algo_solver
    case 'sara'
        for iEffCh = 1:nEffectiveChans
            noise_std = []; noise_full_vect = []; dieloaded = []; Lambda = []; % init to avoid errors
            for iCh = 1:numel(effChans2Image{iEffCh})
                for idSet = 1:nDataSets
                    %% get data & nW (natural weights) & Dies (if available) for now
                    if flag_apply_imaging_weights
                        dataloaded = load(dataFilename(idSet, effChans2Image{iEffCh}(iCh)), 'nW', 'y', 'nWimag');
                        try y{iEffCh, 1}{iCh}{idSet} = double(dataloaded.y(:)) .* double(dataloaded.nW(:)) .* (double(dataloaded.nWimag(:))); % data whitening
                        catch ; fprintf('\nWARNING: imaging weights could not be applied.');
                            flag_apply_imaging_weights = false;
                            y{iEffCh, 1}{iCh}{idSet} = double(dataloaded.y(:)) .* double(dataloaded.nW(:));
                        end
                    else
                        dataloaded = load(dataFilename(idSet, effChans2Image{iEffCh}(iCh)), 'nW', 'y');
                        y{iEffCh, 1}{iCh}{idSet} = double(dataloaded.y(:)) .* double(dataloaded.nW(:));
                    end
                    % data whitening by applying the square root of the weights
                    clear dataloaded;
                    %% get noise STD of the physical channels & l2 bounds (if no visibility gridding)

                    tmp_nmeas = numel(y{iEffCh, 1}{iCh}{idSet});
                    chi_square_bound = sqrt(tmp_nmeas + 2 * sqrt(tmp_nmeas));
                    if ~flag_noise_theoretical

                        tmpfile = param_preproc_files.noise_std(idSet, effChans2Image{iEffCh}(iCh));
                        if isfile(tmpfile)
                            try noiseEffChansloaded = load(tmpfile, 'sigma');
                                if ~flag_visibility_gridding
                                    % case: no visibility gridding:
                                    l2bounds{iEffCh, 1}{iCh}{idSet} =  chi_square_bound * noiseEffChansloaded.sigma;
                                else
                                    % case: visibility gridding: get std of the noise data to be
                                    % used later for the  estimation of the noise std
                                    % after visibility gridding
                                    noise_std{iEffCh, 1}{iCh, 1}{idSet, 1} = noiseEffChansloaded.sigma;
                                end
                                noise_full_vect = [noise_full_vect; noiseEffChansloaded.sigma * (randn(tmp_nmeas, 1) + 1i * randn(tmp_nmeas, 1)) / sqrt(2)];
                            catch
                                flag_noise_theoretical = 1;
                                l2bounds{iEffCh, 1}{iCh}{idSet} = chi_square_bound;
                                fprintf('\nWARNING: failed to load `sigma` from the noise file.');
                            end
                        else; fprintf('\nWARNING: noise statistics file not found: %s ', tmpfile);
                            flag_noise_theoretical = 1;
                        end
                    elseif  ~flag_visibility_gridding
                        % case: no visibility gridding: compute l2bounds if
                        % not available
                        l2bounds{iEffCh, 1}{iCh}{idSet} = chi_square_bound;
                    end

                end
            end
            %% noise STD of the effective channel
            if flag_noise_theoretical
                global_sigma_noise(iEffCh, 1) = 1;
                fprintf('\nINFO: theoretical noise assumed: sigma %f.\n', global_sigma_noise(iEffCh, 1));
            else
                global_sigma_noise(iEffCh, 1) = std(noise_full_vect);
                clear noise_full_vect;
                fprintf('\nINFO: user input: corrected sigma %f.\n', global_sigma_noise(iEffCh, 1));
            end
        end
        %% get measurement op.
        if flag_visibility_gridding
            % one data block & one channel after visibility gridding (DR)
            % ! `G` is the holographic matrix
            % get the reduced data vector `y`, the weighting
            % matrix `Lambda` involved in visibility gridding and
            % the l2 bounds computed from a realisation of a random
            % variable (variance from user or theoretical assumed)
            [A, At, G, W, aW, Lambda, y, noise] = util_load_data_gen_dr_measurement_operator(y, dataFilename, effChans2Image, ...
                nDataSets, Nx, Ny, param_nufft, param_wproj, noise_std, param_preproc_files);
            % get the l2-bounds after gridding the noise
            for iEffCh = 1:nEffectiveChans
                l2bounds{iEffCh}{1} = noise.l2bounds{iEffCh}{1};
                global_sigma_noise(iEffCh, 1) = noise.sigma{iEffCh};
            end
        else
            % re-structure data & l2 bounds; collapse cells
            for iEffCh = 1:nEffectiveChans
                y{iEffCh, 1} = vertcat(y{iEffCh, 1}{:}); y{iEffCh, 1} = y{iEffCh, 1}(:);
                l2bounds{iEffCh, 1} = vertcat(l2bounds{iEffCh, 1}{:}); l2bounds{iEffCh, 1} = l2bounds{iEffCh, 1}(:);
            end
            [A, At, G, W, aW] = util_load_data_gen_measurement_operator(dataFilename, effChans2Image, nDataSets, ...
                Nx, Ny, param_nufft, param_wproj, param_precond, param_preproc_files.dde);
        end
        % get dirty image
        dirty = compute_residual_images(zeros(Ny, Nx), y, A, At, G, W, flag_visibility_gridding, Lambda);
        % setup parpool
        try delete(gcp('nocreate')); end
        hpc_cluster = util_set_parpool('sara', ncores_data, 1, parcluster);

    otherwise % 'hs' or 'fhs'
        %% setup parpool
        try delete(gcp('nocreate')); end
        hpc_cluster = util_set_parpool(algo_solver, ncores_data, Qx * Qy, parcluster);

        % create the measurement operator operator in parallel
        data_first_labindex_minus_1 = Qx * Qy * strcmp(algo_solver, 'fhs'); % index of last prior worker
        spmd
            Lambda = [];
            if labindex > data_first_labindex_minus_1 % define operator on data workers only
                local_fc = (freqRangeCores(labindex - data_first_labindex_minus_1, 1):freqRangeCores(labindex - data_first_labindex_minus_1, 2));
                effChans2Image_lab = (effChans2Image(local_fc));
                for ifc  = 1:numel(local_fc)
                    noise_full_vect = []; dieloaded = []; noise_std = [];
                    for iCh = 1:numel(effChans2Image_lab{ifc})
                        fprintf('\nChannel %d,', effChans2Image_lab{ifc}(iCh));
                        for idSet = 1:nDataSets
                            % get data & apply nW (natural weights)
                            if flag_apply_imaging_weights
                                dataloaded = load(dataFilename(idSet, effChans2Image_lab{ifc}(iCh)), 'nW', 'y', 'frequency', 'nWimag');
                                try y{ifc, 1}{iCh}{idSet} = double(dataloaded.y(:)) .*  double(dataloaded.nW(:)) .*  (double(dataloaded.nWimag(:))); % data whitening
                                catch ; fprintf('\nWARNING: imaging weights could not be applied.');
                                    y{ifc, 1}{iCh}{idSet} = double(dataloaded.y(:)) .*  double(dataloaded.nW(:));
                                end
                            else
                                dataloaded = load(dataFilename(idSet, effChans2Image_lab{ifc}(iCh)), 'nW', 'y', 'frequency');
                                y{ifc, 1}{iCh}{idSet} = double(dataloaded.y(:)) .*  double(dataloaded.nW(:)); % data whitening
                            end

                            frequency =  dataloaded.frequency;

                            dataloaded = [];
                            % get noise STD of the physical channels & l2 bounds (if no visibility gridding)
                            tmp_nmeas = numel(y{ifc, 1}{iCh}{idSet});
                            chi_square_bound = sqrt(tmp_nmeas + 2 * sqrt(tmp_nmeas));
                            if ~flag_noise_theoretical
                                tmpfile = param_preproc_files.noise_std(idSet, effChans2Image_lab{ifc}(iCh));
                                if isfile(tmpfile)
                                    try noiseEffChansloaded = load(tmpfile, 'sigma');
                                        if ~flag_visibility_gridding
                                            % case: no visibility gridding:
                                            l2bounds{ifc, 1}{iCh}{idSet} =  chi_square_bound * noiseEffChansloaded.sigma;
                                        else
                                            % case:  visibility gridding: get std of the noise data to be
                                            % used later for the  estimation of the noise level
                                            % after visibility gridding
                                            noise_std{ifc, 1}{iCh, 1}{idSet, 1} = noiseEffChansloaded.sigma;
                                        end
                                        noise_full_vect = [noise_full_vect; ...
                                            noiseEffChansloaded.sigma * (randn(tmp_nmeas, 1) + 1i * randn(tmp_nmeas, 1)) / sqrt(2)];
                                    catch
                                        flag_noise_theoretical = 1;
                                        fprintf('\nWARNING: failed to load `sigma` from the noise file.');
                                        l2bounds{ifc, 1}{iCh}{idSet} = chi_square_bound;
                                    end
                                else;  fprintf('\nWARNING: noise statistics file not found: %s ', tmpfile);
                                    flag_noise_theoretical = 1;
                                end
                            elseif  ~flag_visibility_gridding
                                % case: no visibility gridding, compute l2bounds if not found
                                l2bounds{ifc, 1}{iCh}{idSet} =  chi_square_bound;
                            end

                        end
                    end
                    % noise STD of the effective channel
                    if flag_noise_theoretical
                        global_sigma_noise_cmpst(ifc, 1) = 1;
                        fprintf('\nINFO: theoretical noise assumed: sigma %f', global_sigma_noise_cmpst(ifc, 1));
                    else
                        global_sigma_noise_cmpst(ifc, 1) = std(noise_full_vect); noise_full_vect = [];
                        fprintf('\nINFO: user input: corrected sigma %f', global_sigma_noise_cmpst(ifc, 1));
                    end

                end
                % build the measurement operator
                if flag_visibility_gridding
                    % one data block & one channel after visibility gridding
                    % ! `G` is the holographic matrix
                    % get the reduced data vector `y`, the weighting
                    % matrix `Lambda` involved in visibility gridding
                    % and the computed l2 bounds from a realisation of a random
                    % variable (std from user or theoretical assumed)
                    fprintf('\nCompute the holographic matrix  ..');
                    [A, At, G, W, aW, Lambda, y, noise] = util_load_data_gen_dr_measurement_operator(y, dataFilename, effChans2Image_lab, ...
                        nDataSets, Nx, Ny, param_nufft, param_wproj, noise_std, param_preproc_files);
                    % get the l2-bounds from the gridded noise vector
                    for ifc = 1:numel(local_fc)
                        l2bounds{ifc}{1} = noise.l2bounds{ifc}{1};
                        global_sigma_noise_cmpst(ifc, 1) = noise.sigma{ifc};
                    end
                    noise_std = [];
                else
                    Lambda = [];
                    % re-structure data; collapse cells
                    for ifc  = 1:numel(local_fc)
                        y{ifc, 1} = vertcat(y{ifc, 1}{:}); y{ifc, 1} =  y{ifc, 1}(:);
                        l2bounds{ifc, 1} = vertcat(l2bounds{ifc, 1}{:}); l2bounds{ifc, 1} = l2bounds{ifc, 1}(:);
                    end
                    % ! `G` is the degridding matrix
                    [A, At, G, W, aW] = util_load_data_gen_measurement_operator(dataFilename, effChans2Image_lab, nDataSets, ...
                        Nx, Ny, param_nufft, param_wproj,  param_precond, param_preproc_files.dde);
                end

                dirty = compute_residual_images(squeeze(zeros(Ny, Nx, numel(y))), y, A, At, ...
                    G, W, flag_visibility_gridding, Lambda);
            end
        end
end
% Free memory
clear  param_wproj param_preproc  local_fc preproc_*;
fprintf('\nData loaded, and measurement operator built, successfully.');
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% save dirty images
try
    subcube_channels = 1:nEffectiveChans;
    switch algo_solver
        case 'sara'
            dirty_file = fullfile(results_path, ...
                strcat('dirty_', algo_solver, ...
                '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
                '_chs', num2str(effChans2Image{1}(1)), '-', num2str(effChans2Image{1}(end)), '.fits'));
            fitswrite(dirty, dirty_file);
        otherwise
            for k = 1:ncores_data
                dirty_file = fullfile(results_path, strcat('dirty_', algo_solver, ...
                    '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
                    '_ch', num2str(k), '.fits'));
                dirty_tmp = dirty{data_worker_id(k)};
                fitswrite(dirty_tmp, dirty_file);
            end
    end
    fprintf('\nDirty images saved, successfully.');
catch
    fprintf('\nWARNING: error encountered when saving dirty images.');
end
clear dirty dirty_tmp;
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% generate l2 bounds
switch algo_solver
    case 'sara'
        % ! to be verified
        % all the variables are stored on the main process for sara
        epsilons = l2bounds(1);
    otherwise
        epsilons = Composite();
        for k = 1:ncores_data
            epsilons{data_worker_id(k)} = l2bounds{data_worker_id(k)};
            if k == 1; global_sigma_noise = global_sigma_noise_cmpst{data_worker_id(k)};
            else; global_sigma_noise = [global_sigma_noise; global_sigma_noise_cmpst{data_worker_id(k)}];
            end
        end; clear global_sigma_noise_cmpst;
end
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% Compute operator norm
subcube_channels = 1:nEffectiveChans;
if strcmp(algo_solver, 'sara')
    if flag_computeOperatorNorm
        fprintf('\nINFO: computing operator''s spectral norm .. ');
        tic;
        [Anorm, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-6, 200, flag_visibility_gridding, Lambda); % AD: changed tolerance to 1e-6 instead of 1e-8
        save(fullfile(results_path, ...
            strcat('Anorm_', algo_solver, ...
            '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
            '_chs', num2str(effChans2Image{1}(1)), '-', num2str(effChans2Image{1}(end)), '.mat')), ...
            '-v7.3', 'Anorm', 'squared_operator_norm', 'rel_var', ...
            'squared_operator_norm_precond', 'rel_var_precond');
        clear rel_var;
        fprintf('done in %d sec.\n', floor(toc));
    else
        load(fullfile(results_path, ...
            strcat('Anorm_', algo_solver, ...
            '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
            '_chs', num2str(effChans2Image{1}(1)), '-', num2str(effChans2Image{1}(end)), '.mat')), ...
            'Anorm', 'squared_operator_norm_precond', 'squared_operator_norm');
    end
else
    if flag_computeOperatorNorm
        tic;
        fprintf('\nINFO: computing operator''s spectral norm .. ');
        spmd
            if labindex > data_first_labindex_minus_1
                [An, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-8, 200, flag_visibility_gridding, Lambda);
            end
        end

        if flag_saveOperatorNorm
            % save operator norm from the different subcubes into a single .mat
            % file
            fprintf('\nINFO: measurement operator''s norm is saved. \n');
            opnormfile = matfile(fullfile(results_path, strcat('Anorm', ...
                '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
                '_L', num2str(nEffectiveChans), '.mat')), 'Writable', true);
        end
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

        fprintf('done in %d sec.\n', floor(toc));
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
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
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
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% Regularization  and solver parameters
% get the noise level (set regularization parameters to the same value)
% compute sig and sig_bar (get  the "noise level" in "SVD" and
% SARA space) involved in the reweighting scheme

switch algo_solver
    case 'sara'
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
    otherwise
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

% additional solver params.
param_solver.squared_operator_norm = squared_operator_norm;
param_solver.reweighting_sig = sig; % estimate of the noise level in SARA space
if ~strcmp(algo_solver, 'sara')
    param_solver.reweighting_sig_bar = sig_bar; % estimate of the noise level in "SVD" spaces
end
param_solver.nu0 = 1; % bound on the norm of the Identity operator
param_solver.nu1 = 1; % bound on the norm of the operator Psi
param_solver.nu2 = squared_operator_norm_precond; % upper bound on the norm of the measurement operator
if ~strcmp(algo_solver, 'sara')
    param_solver.gamma0 = mu_bar; % regularization parameter nuclear norm
end
param_solver.reweighting_flag = flag_reweighting;
if ~flag_reweighting
    param_solver.reweighting_max_iter = 0;
end
param_solver.gamma = mu; % regularization parameter l21-norm (soft th parameter) ! for SARA, take the value given as an input to the solver
param_solver.cube_id = subcubeInd; % id of the cube to be reconstructed
param_solver.backup_frequency = 1; % PA: checkpoint frequency!

param_solver.alph = gam;
param_solver.alph_bar = gam_bar;
param_solver.adjust_flag_noise = adjust_flag_noise;
% temp filenames
name_checkpoint = fullfile(auxiliary_path, temp_results_name(nEffectiveChans));
name_warmstart = fullfile(auxiliary_path, warm_start(nEffectiveChans));
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% image cube initialisation: load if available otherwise set to 0
switch algo_solver
    case 'sara'
        xinit = zeros(Ny, Nx);
        if ~isempty(filename_init_model)
            fprintf('\nLoading  available image estimate  ...');
            try xinit = fitsread(filename_init_model(effChans2Image{1}(1), effChans2Image{1}(end)));
                if  ~(size(xinit, 1) == Ny && size(xinit, 2) == Nx)
                    xinit = zeros(Ny, Nx);
                    fprintf('\nWARNING: init. image with different dimensions found.');
                else; fprintf('\nINFO: init. image loaded successfully.');
                end
            catch; fprintf('\nWARNING: init. image not found.');
            end
        end
    otherwise
        xinit = zeros(ImageCubeDims);
        if ~isempty(filename_init_model)
            fprintf('\nLoading  available image estimates  ...');
            for  iEffCh =  1:nEffectiveChans
                try  xinit(:, :, iEffCh) = reshape(fitsread(filename_init_model(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end))), N, 1);
                    fprintf('\nINFO: Effective channel ID %d: init. image loaded successfully.', iEffCh);
                catch; fprintf('\nWARNING: Effective channel ID %d: init. image not found or wrong dimensions.', iEffCh);
                end
            end
        end
end
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% Run imaging
if flag_solveMinimization
    if strcmp(algo_solver, 'sara')
        disp('SARA');
        disp('-----------------------------------------');
        % ! in this case, ncores_data corresponds
        % ! to the number of workers for the wavelet transform (9 maximum)
        xsol = sara(y, epsilons, A, At, aW, G, W, Psi, Psit, ...
            param_solver, name_warmstart, name_checkpoint, gam, ...
            flag_visibility_gridding, Lambda, xinit);
        fitswrite(xsol, fullfile(auxiliary_path, strcat('MODEL_', srcName, '_', algo_solver, ...
            '_', num2str(pixelSize), 'asec', ...
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
        switch algo_solver
            case 'hs'
                disp('HyperSARA');
                disp('-----------------------------------------');
                xsol = hyperSARA(y, epsilons, ...
                    A, At, aW, G, W, param_solver, ...
                    ncores_data, dict.basis, dict.nlevel, channel_chunks, ...
                    nEffectiveChans, Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flag_visibility_gridding, Lambda, ...
                    xinit);
            case 'fhs'
                disp('Faceted HyperSARA');
                disp('-----------------------------------------');
                xsol = facetHyperSARA(y, epsilons, ...
                    A, At, aW, G, W, param_solver, Qx, Qy, ncores_data, ...
                    dict.basis, dict.filter_length, dict.nlevel, window_type, ...
                    channel_chunks, nEffectiveChans, overlap_size, ...
                    Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flag_visibility_gridding, Lambda, ...
                    xinit);
            otherwise; error('Unknown solver version.');
        end
        fitswrite(xsol, fullfile(auxiliary_path, strcat('WB_MODEL_', srcName, '_', algo_solver, ...
            '_', num2str(pixelSize), 'asec', ...
            '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
            '_ind', num2str(subcubeInd), ...
            '_gam', num2str(gam), '_gambar', num2str(gam_bar), ...
            '_overlap', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
            '.fits')));
    end
end
