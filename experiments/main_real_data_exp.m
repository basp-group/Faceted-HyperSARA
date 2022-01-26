function main_real_data_exp(image_name, datasetsNames, dataFilename, ...
    subcubeInd, effChans2Image, param_solver, param_global, ...
    param_nufft, param_blocking, param_precond, param_nnls, dict)
% Main script to run the faceted HyperSARA approach on real data.
%
% This script generates synthetic data and runs the SARA, HyperSARA or
% faceted HyperSARA approach to reconstruct an :math:`N \times L` wideband
% image cube.
%
% Parameters
% ----------
% image_name : string
%     Name of the reference synthetic image (from the ``data/`` folder).
% dataSetsNames: cell of string
%     Names of the datasets to be imaged.
% dataFilenames: cell of string
%     Name of the data set files.
% subcube_ind : int
%     Index of the spectral facet to be reconstructed (set to -1 or 0 to
%     deactivate spectral faceting).  AD: is this  still the case ????
% effChans2Image: cell array
%     Indices of the "physical" channels to be concatenated for each 
%     effective channel.
% param_global: struct
%     Global imaging pipeline parameters (see 
%     :mat:func:`experiments.main_input_exp`).
% param_global.algo_version : string (``"sara"``, ``"hs"`` or ``"fhs"``).
%     Selected solver.
% param_global.ncores_data : int
%     Number of cores handlig the data fidelity terms (data cores). For
%     Faceted HyperSARA, the total number of cores used
%     is ``Qx*Qy + ncores_data + 1``. For SARA and HyperSARA, represents
%     the number of cores used for the parallelization.
% param_global.im_pixelSize : double
%     pixel size in arcsec.
% param_global.im_Nx : int
%     Image dimension (x axis, dimension 2).
% param_global.im_Ny : int
%     Image dimension (y axis, dimension 1).
% param_global.facet_Qx : int
%     Number of spatial facets along axis 2 (axis x).
% param_global.facet_Qy : int
%     Number of spatial facets along axis 1 (axis y).
% param_global.facet_Qc : int
%     Number of spectral facets.
% param_global.facet_window_type : string (``"triangular"``, ``"hamming"`` or ``"pc"`` (piecewise-constant))
%     Type of apodization window considered for the faceted nuclear norm
%     prior (FHS solver).
% param_global.facet_overlap_fraction : double[2]
%     Fraction of the total size of a facet overlapping with a neighbour 
%     facet along each axis (y and x) for the faceted low-rankness prior.
%     Will be reset automatically to ``[0, 0]`` if 
%     ``param_global.algo_version = "sara"``
%     or ``"hs"``. Besides, each entry of
%     ``param_global.facet_overlap_fraction`` is reset to 0 the
%     number of facet along the corresponding dimension is equal to 1
%     (i.e., ``param_global.facet_Qy = 1`` and/or
%     ``param_global.facet_Qx = 1``).
% param_global.reg_nReweights : int
%     Maximum number of reweighting steps.
% param_global.reg_gam : double
%     Additional multiplicative factor affecting the joint-sparsity
%     regularization term.
% param_global.reg_gam_bar : double
%     Additional multiplicative factor affecting the low-rankness
%     regularization term.
% param_global.reg_flag_reweight : int
%     Flag to activate re-weighting.
% param_global.reg_flag_homotopy : int
%     Flag to activate  homotopy strategy in the re-weighting.
% param_global.algo_flag_computeOperatorNorm : bool
%     Flag triggering the computation of the (preconditioned) operator norm.
% param_global.algo_flag_solveMinimization : bool
%     Flag triggering the solver (SARA, HS or FHS).
% param_global.algo_flag_dataReduction : bool
%     Flag to activate Data reduction features in the definition of the
%     measurement operator.
% cirrus:
%     hardware: Specify whether the solver runs on cirrus or not (for the
%     creation of the parpool).
%
% Note
% ----
% - DR features still need to be implemented in the main script.
% - Example:
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


% ------------------------------------------------------------------------%
% ------------------------------------------------------------------------%
%% constant
speed_of_light = 299792458;

%% check input params : @PA in light of the new structure of params, many of these fields will be deleted I guess ..
% Image resolution & dimensions
% ! parameters and default values (keep only what is required)
if ~isfield(param_global, 'im_Nx');  param_global.im_Nx = 2048; end
if ~isfield(param_global, 'im_Ny');  param_global.im_Ny = 2048; end
if ~isfield(param_global, 'im_pixelSize');  param_global.im_pixelSize = []; end


% Prior: wideband facet related
if ~isfield(param_global, 'facet_Qx');  param_global.facet_Qx =  floor(param_global.facet_Nx / 256); end
if ~isfield(param_global, 'facet_Qy');  param_global.facet_Qy =  floor(param_global.facet_Ny / 256); end
if ~isfield(param_global, 'facet_Qc'); param_global.facet_Qc = 1; end
if ~isfield(param_global, 'facet_window_type'); param_global.facet_window_type = 'triangular'; end
if ~isfield(param_global, 'facet_overlap_fraction'); param_global.facet_overlap_fraction = [0.1, 0.1]; end
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
if ~isfield(param_global, 'reg_flag_homotopy'); param_global.reg_flag_homotopy = 0; end

% Data blocks
% if ~isfield(param_global, 'generate_eps_nnls'); param_global.generate_eps_nnls = false; end
% if ~isfield(param_global, 'data_nDataBlk');  param_global.data_nDataBlk = []; end
% if ~isfield(param_global, 'data_sizeDataBlk');  param_global.data_sizeDataBlk = []; end
% if isempty(param_global.data_sizeDataBlk) && isempty (param_global.data_nDataBlk)
%     param_global.data_sizeDataBlk = 2e5;
% end

% Algo
if ~isfield(param_global, 'algo_version'); param_global.algo_version = 'fhs'; end
if ~isfield(param_global, 'algo_flag_solveMinimization'); param_global.algo_flag_solveMinimization = 1; end
if ~isfield(param_global, 'algo_flag_computeOperatorNorm'); param_global.algo_flag_computeOperatorNorm = 1; end
if ~isfield(param_global, 'algo_ncores_data'); param_global.algo_ncores_data = numel(effChans2Image);  end

% Measurement operator
if ~isfield(param_global, 'measop_flag_dataReduction'); param_global.measop_flag_dataReduction = 0; end
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

% Hardware
hpc = param_global.hardware;
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
    reg_nReweights = param_global.reg_nReweights; % @PA: why commented?
    flag_homotopy = param_global.reg_flag_homotopy;
end
% % Data blocks
% nDataBlk = param_global.data_nDataBlk;
% szDataBlk = param_global.data_sizeDataBlk;
% Meas. op.
flag_dataReduction = param_global.measop_flag_dataReduction; % data reduction
param_wproj.do = param_global.measop_flag_wproj; % w projection
if isempty(param_global.measop_flag_wproj)
    param_wproj.do =0;
end
param_wproj.CEnergyL2 = param_global.measop_wprojCEnergyL2; % w projection
param_wproj.GEnergyL2 = param_global.measop_wprojGEnergyL2; % w projection
% Algo
algo_version = param_global.algo_version;
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

% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% Paths
format compact;
project_dir = param_global.main_dir;
fprintf('\nMain project dir. is %s: ', project_dir);
current_dir = pwd;
if strcmp(project_dir, current_dir)
    current_dir = [project_dir, 'experiments', filesep];
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
addpath([current_dir, filesep, 'real_data']);% AD: will be ommited later
%
if strcmp(algo_version, 'sara')
    addpath([project_dir, filesep, 'src', filesep, 'sara']);
elseif strcmp(algo_version, 'hs')
    addpath([project_dir, filesep, 'src', filesep, 'hs', filesep]);
elseif strcmp(algo_version, 'fhs')
    addpath([project_dir, filesep, 'src', filesep, 'fhs', filesep]);
end
% setting paths to results and reference image cube
results_path = fullfile('results', image_name);
mkdir(results_path);
auxiliary_path = fullfile(results_path, algo_version);
mkdir(auxiliary_path);
fprintf('\nINFO: results will be saved in: \n %s .\n', auxiliary_path);
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% info
disp('Imaging configuration');
disp(['Source name: ', image_name]);
disp(['Image cube size: ', num2str(Ny), ' x ', num2str(Nx), ' x ', num2str(numel(effChans2Image))]);
disp(['Algorithm version: ', algo_version]);
disp(['Number of facets Qy x Qx : ', num2str(Qy), ' x ', num2str(Qx)]);
disp(['Number of spectral facets Qc : ', num2str(Qc)]);
if ~strcmp(algo_version, 'sara')
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
if strcmp(imResolution,'nominal')
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
    pixelSize = (180 / pi) * 3600 / (2 * spatialBandwidth);
    fprintf('\nINFO: the pixelsize is set to %f arcsec, that is 2x nominal resolution at the highest freq.', pixelSize);
end
halfSpatialBandwidth = (180 / pi) * 3600 / (pixelSize) / 2;
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% image cube initialisation: load if available otherwise set to 0
switch algo_version
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
            catch; fprintf('\nWARNING: init. image not found.');
            end
        end
    otherwise
        xinit = zeros(N, nEffectiveChans);
        if ~isempty(param_preproc.filename_model)
            fprintf('\nLoading  available image estimate  ...');
            for  iEffCh =  1:nEffectiveChans
                try  xinit(:, iEffCh) = reshape(fitsread(param_preproc.filename_model(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end))), N, 1);
                    fprintf('\nINFO: init. image found.');
                catch; fprintf('\nWARNING: init. image not found.');
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
        l2EffChansloaded = load(filename_l2bounds(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end)), 'sigmac', 'l2bounds');
        if ~flag_dataReduction
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
            catch; fprintf('\nWARNING: l2 bounds not found. \nChi squared distribution is assumed to determine the l2 bounds.\n');
                flag_l2bounds_compute = 1;
            end
        end
        if ~flag_l2bounds_compute
            try % get noise std
                if ~flag_dataReduction
                    global_sigma_noise(iEffCh, 1) = full(l2EffChansloaded.sigmac);
                else;  global_sigma_noise(iEffCh, 1) = full(l2EffChansloaded.sigmac.gridded);
                end
            catch;  global_sigma_noise(iEffCh, 1) = 1;
                flag_l2bounds_compute = 1; % assuming a chi squared dist.
                fprintf('\nWARNING: global sigma noise not found, will assume white data.\n');
            end
        end
    end
end
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% Auxiliary function needed to select the appropriate workers
% (only needed for 'hs' and 'fhs' algorithms)
switch algo_version
    case 'sara'
        data_worker_id = @(k) k;
    case 'hs'
        data_worker_id = @(k) k;
    case 'fhs'
        data_worker_id = @(k) k + Qx * Qy;
    otherwise
        error('Undefined algo_version (available options: ''sara'',''hs'',''fhs'')');
end
%% Get faceting parameter (spectral + spatial)
% fix faceting parameters in case these are not consistent with the
% selected algorithm
if strcmp(algo_version, 'sara')
    window_type = 'none';
    Qc = 1; % nEffectiveChans;
    Qx = 1;
    Qy = 1;
elseif strcmp(algo_version, 'hs')
    window_type = 'none';
    Qx = 1;
    Qy = 1;
end
Q = Qx * Qy;

% convert fraction of overlap between consecutive facets into a number of pixels
overlap_size = get_overlap_size([Ny, Nx], [Qy, Qx], overlap_fraction);
if ~strcmp(algo_version, 'sara')
    fprintf('\nINFO: faceting: number of pixels in overlap: [%s]\n', strjoin(strsplit(num2str(overlap_size))));
end
% index of channels from the subcube to be handled on each data worker
if strcmp(algo_version, 'sara')
    ncores = numel(dict.basis );
    ncores_data = ncores;
end
freqRangeCores = split_range(ncores_data, nEffectiveChans);
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
delete(gcp('nocreate'));
if strcmp(algo_version, 'sara')
    hpc = 'local';
end
cirrus_cluster = util_set_parpool_dev(algo_version, ncores_data, Qx * Qy, strcmp(hpc, 'cirrus'));
% -------------------------------------------------------------------------%
% -------------------------------------------------------------------------%
%% Setup measurement operator and load data
switch algo_version
    case 'sara'
        for iEffCh = 1:nEffectiveChans
            
            %% if die/dde calib, get dies/ddes
            if flag_calib.die
                dieloaded = load(param_preproc.filename_die(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end)), 'DIEs');
            elseif flag_calib.dde
                ddeloaded = load(param_preproc.filename_dde(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end)), 'DDEs');
            end
            %% if data reduction and residual data available, load them
            if flag_dataReduction
                DR_residualdataloaded = load(filename_l2bounds(effChans2Image{iEffCh}(1), effChans2Image{iEffCh}(end)), 'RESIDUAL_DATA');
            end
            for iCh = 1:numel(effChans2Image{iEffCh})
                for idSet = 1:nDataSets
                    dataloaded = load(dataFilename(idSet, effChans2Image{iEffCh}(iCh)), 'u', 'v', 'w', 'nW', 'y','frequency');
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
                    
                    if flag_calib.die
                        % if die calib, correct data
                        y{iEffCh, 1}{iCh}{idSet} = y{iEffCh, 1}{iCh}{idSet} ./ dieloaded.DIEs{idSet, iCh};
                        dieloaded.DIEs{idSet, iCh} = [];
                    elseif flag_calib.dde
                        % dde calib
                        ddes{iEffCh, 1}{iCh}{idSet} = ddeloaded.DDEs{idSet, iCh};
                        ddeloaded.DDEs{idSet, iCh} = [];
                    end
                    if flag_dataReduction
                        % if residual data are available, use them to get
                        % l2 bounds.
                        try preproc_dr_residuals{iEffCh, 1}{iCh}{idSet} = DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh};
                            DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh} = [];
                        catch,  preproc_dr_residuals = [];
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
            else, ddes = [];
            end
            
            if flag_dataReduction
                % if residual data are available
                if exist('preproc_dr_residuals', 'var')
                    preproc_dr_residuals{iEffCh, 1} = vertcat(preproc_dr_residuals{iEffCh, 1}{:});
                    preproc_dr_residuals{iEffCh, 1} = preproc_dr_residuals{iEffCh, 1}(:);
                end
            end
        end
        
        if flag_dataReduction
            % (DR) one data block & one channel 
            [A, At, G, W, aW, Sigma, y, noise] = util_gen_dr_measurement_operator_dev(y, u, v, w, nW, ...
                1, Nx, Ny, param_nufft, param_wproj,param_precond, param_blocking,preproc_dr_residuals);
            % (DR) get the computed l2 bounds from the continuous residual
            % data if available
            if isfield(noise,'l2bounds') && isfield(noise,'sigma')
                for iEffCh = 1:nEffectiveChans
                    l2bounds{iEffCh}{1} = noise.l2bounds{iEffCh}{1};
                    global_sigma_noise(iEffCh, 1) = noise.sigma{iEffCh};
                end
            elseif ~isempty(preproc_dr_residuals)
                fprintf('\nWARNING: computing l2 bounds from the given residual data failed. ');
            end
        else
            [A, At, G, W, aW] = util_gen_measurement_operator_dev(u, v, w, nW, ...
                param_precond, param_blocking, 1, Nx, Ny, param_nufft, param_wproj, ddes);
            Sigma = [];
        end; clear u v w nW;
        
    otherwise % 'hs' or 'fhs'
        % create the measurement operator operator in parallel (depending on
        % the algorithm used)
        Sigma = Composite();
        if strcmp(algo_version, 'hs')
            spmd
                ddes = []; preproc_dr_residuals = [];
                local_fc = (freqRangeCores(labindex, 1):freqRangeCores(labindex, 2));
                effChans2Image_lab = (effChans2Image(local_fc));
                for ifc  = 1:numel(local_fc)
                    %% if die/ddes calib, get dies/ddes
                    if flag_calib.die
                        dieloaded = load(param_preproc.filename_die(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'DIEs');
                    elseif flag_calib.dde
                        ddeloaded = load(param_preproc.filename_dde(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'DDEs');
                    end
                    %% if data reduction and residual data available, load them
                    if flag_dataReduction
                        DR_residualdataloaded = load(filename_l2bounds(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'RESIDUAL_DATA');
                    end
                    for iCh = 1:numel(effChans2Image_lab{ifc})
                        for idSet = 1:nDataSets
                            %% load data
                            dataloaded = load(dataFilename(idSet, effChans2Image_lab{ifc}(iCh)), 'u', 'v', 'w', 'nW', 'y','frequency');
                            % u v  are in units of the wavelength and will be normalised between [-pi,pi] for the NUFFT
                            u{ifc, 1}{iCh}{idSet} = double(dataloaded.u(:)) * pi / halfSpatialBandwidth;
                            dataloaded.u = [];
                            v{ifc, 1}{iCh}{idSet} = -double(dataloaded.v(:)) * pi / halfSpatialBandwidth;
                            dataloaded.v = [];
                            w{ifc, 1}{iCh}{idSet} = double(dataloaded.w(:));
                            dataloaded.w = [];
                            nW{ifc, 1}{iCh}{idSet} = double(dataloaded.nW(:)); % sqrt(natural weights)
                            dataloaded.nW = [];
                            y{ifc, 1}{iCh}{idSet} = double(dataloaded.y(:)) .* nW{ifc, 1}{iCh}{idSet};% data whitening
                            dataloaded.y = [];
                            % Output of the pre-processing step if any
                            if flag_calib.die
                                % if die calib, correct data
                                y{ifc, 1}{iCh}{idSet} = y{ifc, 1}{iCh}{idSet} ./ dieloaded.DIEs{idSet, iCh};
                                dieloaded.DIEs{idSet, iCh} = [];
                            elseif flag_calib.dde
                                % if dde calib, get ddes
                                ddes{ifc, 1}{iCh}{idSet} = ddeloaded.DDEs{idSet, iCh};
                                ddeloaded.DDEs{idSet, iCh} = [];
                            end
                            if flag_dataReduction
                                % if residual data are available, use them to get
                                % l2 bounds.
                                try preproc_dr_residuals{ifc, 1}{iCh}{idSet} = DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh};
                                    DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh} = [];
                                catch, preproc_dr_residuals = [];
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
                    
                    if flag_dataReduction
                        try preproc_dr_residuals{ifc, 1} = vertcat(preproc_dr_residuals{ifc, 1}{:});
                            preproc_dr_residuals{ifc, 1} = preproc_dr_residuals{ifc, 1}(:);
                        catch; preproc_dr_residuals = [];
                        end
                    end
                end
                
                if flag_dataReduction
                    % ! define Sigma (weight matrix involved in DR)
                    % ! define G as the holographic matrix
                    fprintf('\nCompute the holographic matrix H .. \n');
                    [A, At, G, W, aW, Sigma, y] = util_gen_dr_measurement_operator_dev(y, u, v, w, nW, ...
                         numel(local_fc), Nx, Ny, param_nufft, param_wproj, param_precond, param_blocking, preproc_dr_residuals, ddes);
                else
                    % ! ideally, simplify irt nufft interface to do so
                    [A, At, G, W, aW] = util_gen_measurement_operator_dev(u, v, w, nW, ...
                         numel(local_fc), Nx, Ny, param_nufft, param_wproj, param_precond, param_blocking , ddes);
                    Sigma = [];
                end
                u = []; v = []; w = []; nW = [];
                dirac = sparse(Ny * 0.5 + 1, Nx * 0.5 + 1, 1, Ny, Nx);
                for l = 1:numel(G)
                    F = HS_forward_operator_G(full(dirac), G(l), W(l), A, flag_dataReduction, Sigma);
                    psf_peak = max(max(HS_adjoint_operator_G(F, G(l), W(l), At, Ny, Nx, flag_dataReduction, Sigma)));
                    fprintf('\nLab %d: INFO: peak of PSF (for residual normalisation): %f', labindex, psf_peak);
                end
            end
        else
            Sigma = Composite();
            spmd
                ddes = []; preproc_dr_residuals = [];
                % define operator on data workers only
                if labindex > Q
                    local_fc = (freqRangeCores(labindex - Q, 1):freqRangeCores(labindex - Q, 2));
                    effChans2Image_lab = (effChans2Image(local_fc));
                    for ifc  = 1:numel(local_fc)
                        %% if die/ddes calib, get dies/ddes
                        if flag_calib.die
                            dieloaded = load(param_preproc.filename_die(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'DIEs');
                        elseif flag_calib.dde
                            ddeloaded = load(param_preproc.filename_dde(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'DDEs');
                        end
                        %% if data reduction and residual data available, load them
                        if flag_dataReduction
                            DR_residualdataloaded = load(filename_l2bounds(effChans2Image{ifc}(1), effChans2Image{ifc}(end)), 'RESIDUAL_DATA');
                        end
                        for iCh = 1:numel(effChans2Image_lab{ifc})
                            for idSet = 1:nDataSets
                                %% load data
                                dataloaded = load(dataFilename(idSet, effChans2Image_lab{ifc}(iCh)), 'u', 'v', 'w', 'nW', 'y','frequency');
                                % u v  are in units of the wavelength and will be normalised between [-pi,pi] for the NUFFT
                                u{ifc, 1}{iCh}{idSet} = double(dataloaded.u(:)) * pi / halfSpatialBandwidth;
                                dataloaded.u = [];
                                v{ifc, 1}{iCh}{idSet} = -double(dataloaded.v(:)) * pi / halfSpatialBandwidth;
                                dataloaded.v = [];
                                w{ifc, 1}{iCh}{idSet} = double(dataloaded.w(:));
                                dataloaded.w = [];
                                nW{ifc, 1}{iCh}{idSet} = double(dataloaded.nW(:)); % sqrt(natural weights)
                                dataloaded.nW = [];
                                y{ifc, 1}{iCh}{idSet} = double(dataloaded.y(:)) .* nW{ifc, 1}{iCh}{idSet};% data whitening
                                dataloaded.y = [];
                                % Output of the pre-processing step if any
                                if flag_calib.die
                                    % if die calib, correct data
                                    y{ifc, 1}{iCh}{idSet} = y{ifc, 1}{iCh}{idSet} ./ dieloaded.DIEs{idSet, iCh};
                                    dieloaded.DIEs{idSet, iCh} = [];
                                elseif flag_calib.dde
                                    % if dde calib, get ddes
                                    ddes{ifc, 1}{iCh}{idSet} = ddeloaded.DDEs{idSet, iCh};
                                    ddeloaded.DDEs{idSet, iCh} = [];
                                end
                                if flag_dataReduction
                                    % if residual data are available, use them to get
                                    % l2 bounds.
                                    try preproc_dr_residuals{ifc, 1}{iCh}{idSet} = DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh};
                                        DR_residualdataloaded.RESIDUAL_DATA{idSet, iCh} = [];
                                    catch, preproc_dr_residuals = [];
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
                        
                        if flag_dataReduction
                            try preproc_dr_residuals{ifc, 1} = vertcat(preproc_dr_residuals{ifc, 1}{:});
                                preproc_dr_residuals{ifc, 1} = preproc_dr_residuals{ifc, 1}(:);
                            catch; preproc_dr_residuals = [];
                            end
                        end
                    end
                    
                    if flag_dataReduction
                        % ! define Sigma (weight matrix involved in DR)
                        % ! define G as the holographic matrix
                        fprintf('\nCompute the holographic matrix H ..\n');
                        [A, At, G, W, aW, Sigma, y] = util_gen_dr_measurement_operator_dev(y, u, v, w, nW, ...
                            numel(local_fc), Nx, Ny, param_nufft, param_wproj, param_precond, param_blocking, preproc_dr_residuals, ddes);
                    else
                        % ! ideally, simplify irt nufft interface to do so
                        [A, At, G, W, aW] = util_gen_measurement_operator_dev(u, v, w, nW, ...
                            numel(local_fc), Nx, Ny, param_nufft, param_wproj,  param_precond, param_blocking, ddes);
                        Sigma = [];
                    end; u = []; v = []; w = []; nW = [];
                    
                    dirac = sparse(Ny * 0.5 + 1, Nx * 0.5 + 1, 1, Ny, Nx);
                    for l = 1:numel(G)
                        F = HS_forward_operator_G(full(dirac), G(l), W(l), A, flag_dataReduction, Sigma);
                        psf_peak = max(max(HS_adjoint_operator_G(F, G(l), W(l), At, Ny, Nx, flag_dataReduction, Sigma)));
                        fprintf('\nLab %d: INFO: peak of PSF (for residual normalisation): %f', labindex, psf_peak);
                    end
                    
                end
            end
        end; clear local_fc  u v w nW dataSpWinloaded;
        
end

clear  resyCmpst param_wproj param_preproc  ; %% Free memory

%% load l2 bounds (generate only full spectral dataset)
% only generatr data in 'hs' or 'fhs' configuration (otherwise, load the data)
% datafile = matfile(fullfile(results_path,filename_l2bounds));
fprintf('\nData loaded successfully.');
switch algo_version
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
if strcmp(algo_version, 'sara')
    if flag_computeOperatorNorm
        [Anorm, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-6, 200, flag_dataReduction, Sigma); % AD: changed tolerance to 1e-6 instead of 1e-8
        save(fullfile(results_path, ...
            strcat('Anorm_', algo_version, ...
            '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
            '_chs', num2str(effChans2Image{1}(1)), '-', num2str(effChans2Image{1}(end)), '.mat')), ...
            '-v7.3', 'Anorm', 'squared_operator_norm', 'rel_var', ...
            'squared_operator_norm_precond', 'rel_var_precond');
        clear rel_var;
    else
        load(fullfile(results_path, ...
            strcat('Anorm_', algo_version, ...
            '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
            '_chs', num2str(effChans2Image{1}(1)), '-', num2str(effChans2Image{1}(end)), '.mat')), ...
            'Anorm', 'squared_operator_norm_precond', 'squared_operator_norm');
    end
else
    if flag_computeOperatorNorm
        spmd
            if labindex > Qx * Qy * strcmp(algo_version, 'fhs')
                [An, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-8, 200, flag_dataReduction, Sigma);
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

if strcmp(algo_version, 'sara') % AD
    temp_results_name = @(nEffectiveChans) strcat(image_name, '_', ...
        algo_version, '_', num2str(pixelSize), 'asec', ...
        '_Ny', num2str(Ny), '_Nx', num2str(Nx), ...
        '_ind', num2str(subcubeInd), '_chs', num2str(effChans2Image{1}(1)), '-', num2str(effChans2Image{1}(end)), ...
        '_g', num2str(gam), '_gb', num2str(gam_bar), '_');
else
    temp_results_name = @(nEffectiveChans) strcat(image_name, '_', ...
        algo_version, '_', num2str(pixelSize), 'asec', ...
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

if strcmp(algo_version, 'sara')
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
    fprintf('Algo: %s, alpha = %.4e, mu = %.4e, sig = %.4e\n', algo_version, gam, mu, sig);
end

if strcmp(algo_version, 'hs') || strcmp(algo_version, 'fhs')
    % noise level / regularization parameter
    [sig, sig_bar, mu_chi, sig_chi, sig_sara] = ...
        compute_noise_level(Ny, Nx, nEffectiveChans, global_sigma_noise(:), ...
        algo_version, Qx, Qy, overlap_size, squared_operator_norm(:));
    % apply multiplicative factor for the regularization parameters (if needed)
    mu_bar = gam_bar * sig_bar;
    mu = gam * sig;
    fprintf('mu_chi = %.4e, sig_chi = %.4e, sig_sara = %.4e\n', mu_chi, sig_chi, sig_sara);
    fprintf('Noise levels: sig = %.4e, sig_bar = [%.4e, %.4e]\n', sig, min(sig_bar), max(sig_bar));
    fprintf('Additional multiplicative actors gam = %.4e, gam_bar = %.4e\n', gam, gam_bar);
    fprintf('Regularization parameters: mu = %.4e, mu_bar = %.4e\n', mu, mu_bar);
    fprintf('Algo: %s, gam = %.4e, gam_bar = %.4e, mu = %.4e, mu_bar = [%.4e, %.4e]\n', algo_version, gam, gam_bar, mu, min(mu_bar), max(mu_bar));
end

% moved here for transprancy: define parameters for the solver
% List of default solver-specific parameters (reweighting, pdfb, ellispsoid
% prjection, epsilon update scheme).
% * general
% estimate of the noise level in SARA space
param_solver.reweighting_sig = sig;
if ~strcmp(algo_version, 'sara');  param_solver.reweighting_sig_bar = sig_bar; % estimate of the noise level in "SVD" spaces
end
param_solver.nu0 = 1; % bound on the norm of the Identity operator
param_solver.nu1 = 1; % bound on the norm of the operator Psi
param_solver.nu2 = squared_operator_norm_precond; % upper bound on the norm of the measurement operator
if ~strcmp(algo_version, 'sara'); param_solver.gamma0 = mu_bar; % regularization parameter nuclear norm
end
param_solver.gamma = mu; % regularization parameter l21-norm (soft th parameter) ! for SARA, take the value given as an input to the solver
param_solver.cube_id = subcubeInd; % id of the cube to be reconstructed
param_solver.backup_frequency = 1; % PA: checkpoint frequency! AD :????
try param_solver.flag_homotopy = flag_homotopy; % flag homotopy strategy
catch; param_solver.flag_homotopy = 0;
end
param_solver.alph = gam;
param_solver.alph_bar = gam_bar;
% temp filenames
name_checkpoint = fullfile(auxiliary_path, temp_results_name(nEffectiveChans));
name_warmstart = fullfile(auxiliary_path, warm_start(nEffectiveChans));
%% Solver
if flag_solveMinimization
    if strcmp(algo_version, 'sara')
        disp('SARA');
        disp('-----------------------------------------');
        % ! in this case, ncores_data corresponds
        % ! to the number of workers for the wavelet transform (9 maximum)
        xsol = sara(y, epsilons, A, At, aW, G, W, Psi, Psit, ...
            param_solver, name_warmstart, name_checkpoint, gam, ...
            flag_dataReduction, Sigma, xinit);
        fitswrite(xsol, fullfile(auxiliary_path, strcat('x_', image_name, '_', algo_version, ...
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
        switch algo_version
            case 'hs'
                disp('HyperSARA');
                disp('-----------------------------------------');
                xsol = hyperSARA(y, epsilons, ...
                    A, At, aW, G, W, param_solver, ...
                    ncores_data, dict.basis, dict.nlevel, channel_chunks, ...
                    nEffectiveChans, Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flag_dataReduction, Sigma, ...
                    xinit);
            case 'fhs'
                disp('Faceted HyperSARA');
                disp('-----------------------------------------');
                % ---
                % ! test: load ground truth image (for debugging purposes)
                exp_type = param_global.exp_type;
                switch exp_type
                    case "spatial"
                        image_name = 'cygASband_Cube_1024_2048_20';
                        spectral_downsampling = 1;
                        spatial_downsampling = 1;
                    case "spectral"
                        image_name = 'cygASband_Cube_256_512_100';
                        spectral_downsampling = 1;
                        spatial_downsampling = 1;
                    case "test"
                        image_name = 'cygASband_Cube_512_1024_20';
                        spectral_downsampling = 10;
                        spatial_downsampling = 2;
                    case "local_test"
                        image_name = 'cygASband_Cube_256_512_100';
                        spectral_downsampling = 25;
                        spatial_downsampling = 1;
                        coverage_path = "data/vla_7.95h_dt10s.uvw256.mat";
                    case "old_local_test"
                        image_name = 'cubeW28';
                        spectral_downsampling = 20;
                        spatial_downsampling = 4;
                        coverage_path = "data/vla_7.95h_dt10s.uvw256.mat";
                    otherwise
                        error("Unknown experiment type");
                end
                reference_cube_path = fullfile('../data', strcat(image_name, '.fits'));
                info        = fitsinfo(reference_cube_path);
                rowend      = info.PrimaryData.Size(1);
                colend      = info.PrimaryData.Size(2);
                sliceend    = info.PrimaryData.Size(3);
                if strcmp(algo_version, 'sara')
                    x0 = fitsread(reference_cube_path, 'primary', ...
                        'Info', info, ...
                        'PixelRegion', {[1 spatial_downsampling rowend], ...
                        [1 spatial_downsampling colend], ...
                        ind});
                else
                    x0 = fitsread(reference_cube_path, 'primary', ...
                        'Info', info, ...
                        'PixelRegion', {[1 spatial_downsampling rowend], ...
                        [1 spatial_downsampling colend], ...
                        [1 spectral_downsampling sliceend]});
                end
                [Ny, Nx, nchans] = size(x0);
                N = Nx * Ny;
                x0 = reshape(x0, [N, nchans]);
                % ---
                
                xsol = facetHyperSARA(y, epsilons, ...
                    A, At, aW, G, W, param_solver, Qx, Qy, ncores_data, ...
                    dict.basis, dict.filter_length, dict.nlevel, window_type, ...
                    channel_chunks, nEffectiveChans, overlap_size, gam, gam_bar, ...
                    Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flag_dataReduction, Sigma, ...
                    xinit, x0);
            otherwise; error('Unknown solver version.');
        end
        fitswrite(xsol, fullfile(auxiliary_path, strcat('x_', image_name, '_', algo_version, ...
            '_', num2str(pixelSize), 'asec', ...
            '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
            '_ind', num2str(subcubeInd), ...
            '_gam', num2str(gam), '_gambar', num2str(gam_bar), ...
            '_overlap', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
            '.fits')));
    end
end