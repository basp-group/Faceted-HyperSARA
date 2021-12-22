function main_real_data_dev_exp2a(image_name, data_file,subcube_ind,spwins2image, ...
    algo_version, ncores_data, param_global, param_reg,input_flags, cirrus)
% Main script to run the faceted HyperSARA approach on synthetic data.
%
% This script generates synthetic data and runs the SARA, HyperSARA or
% faceted HyperSARA approach to reconstruct an :math:`N \times L` wideband
% image cube.
%
% Parameters
% ----------
% image_name : string
%     Name of the reference synthetic image (from the data/ folder).
% nChannels : int
%     Number of spectral channels considered.
% Qx : int
%     Number of spatial facets along axis 2 (x).
% Qy : int
%     Number of spatial facets along axis 1 (y)
% Qc : int
%     Number of spectral facets.
% algo_version : string ('sara', 'hs' or 'fhs')
%     Selected solver.
% window_type : string ('triangular', 'hamming' or 'pc' (piecewise-constant))
%     Type of apodization window considered for the faceted nuclear norm
%     prior (FHS solver).
% ncores_data : int (AD: this should be determined automatically)
%     Number of cores handlig the data fidelity terms ("data cores").
%     For Faceted HyperSARA, the total number of cores used is Qx*Qy +
%     ncores_data + 1. For SARA and HyperSARA, represents the number of
%     cores used for the parallelization.
% subcube_ind : int
%     Index of the spectral facet to be reconstructed (set to -1 to
%     deactivate spectral faceting).
% overlap_fraction : array (1d)
%     Fraction of the total size of a facet overlapping with a neighbour
%     facet.
% nReweights : int
%     Maximum number of reweighting steps.
% coverage_path : string
%     Path and name of the uv-coverage .fits file (w/o file extension).
% gam : double
%     Additional multiplicative factor affecting the joint-sparsity
%     regularization term.
% gam_bar : double
%     Additional multiplicative factor affecting the low-rankness
%     regularization term.
% rw : int
%     [description]
% exp_type : string ('spatial' or 'spectral')
%     Type of the experiment to be reproduced.
% pixelSize : double
%     Coverage superresolution factor.
% isnr : double (AD: not needed for real data)
%     Input SNR used to generate the synthetic visibilities (value in dB).
% flag_generateVisibilities : bool  (AD: not needed for real data)
%     Flag specifying whether the visibilities need to be generated or
%     loaded from an existing .mat file.
% flag_computeOperatorNorm : bool
%     Flag triggering the computation of the (preconditioned) operator
%     norm.
% flag_solveMinimization : bool
%     Flag triggering the solver (SARA, HS or FHS).
% flag_dr : bool
%     Flag to activate DR features in the definition of the measurement
%     operator.
% flag_cirrus : bool
%     Specify whether the solver runs on cirrus or not (for the creation of
%     the parpool).
% flag_homotopy : bool
%     Activate the homotopy strategy within the solver.
%
% ..note::
%    DR features still need to be implemented in the main script.
%
% PARAMETERS FOR DEBUGGING
kernel = 'minmax:tuned'; % 'kaiser' (for real data), 'minmax:tuned'
%% constant
speed_of_light = 299792458;
kernel = 'minmax:tuned'; % 'kaiser' (for real data), 'minmax:tuned'

%% AD: real data

%:AD:Default params to be modified ?
%  algo related
if ~isfield(param_global,'algo_version' ), param_global.algo_version = 'fhs' ; end
% image resolution & dimensions
if ~isfield(param_global, 'Nx'),  param_global.Nx = 2048; end
if ~isfield(param_global, 'Ny'),  param_global.Ny = 2048; end
% if ~isfield(param_global, 'nChannels'),  param_global.nChannels = 1; end

if ~isfield(param_global, 'pixelSize'),   param_global.pixelSize = []; end %in arcsec
if isempty (param_global.pixelSize),imageResolution = 'nominal';
else, imageResolution = 'user_defined';
end
% wideband facet related
if ~isfield(param_global, 'Qx'),  param_global.Qx =  floor(param_global.Nx/256); end
if ~isfield(param_global, 'Qy'),  param_global.Qy =  floor(param_global.Ny/256); end
if ~isfield(param_global, 'Qc'), param_global.Qc = 1; end
if ~isfield(param_global, 'window_type'), param_global.window_type = 'triangular'; end
if ~isfield(param_global, 'overlap_fraction'), param_global.overlap_fraction = [0.5,0.5]; end
% Prior: sparsity dict.
if ~isfield(param_global,'wavelet_level' ), param_global.wavelet_level = 4;  end
if ~isfield(param_global,'wavelet_basis' ), param_global.wavelet_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};  end


% data specific
% data blocks

if ~isfield(param_global, 'nDataBlk'),  param_global.nDataBlk = []; end
if ~isfield(param_global, 'sizeDataBlk'),  param_global.sizeDataBlk = []; end
if isempty(param_global.sizeDataBlk) && isempty (param_global.nDataBlk)
    param_global.sizeDataBlk = 2e5;
end
% w projection
if ~isfield(param_global, 'CEnergyL2') , param_global.CEnergyL2 =1-1e-4; end
if ~isfield(param_global, 'GEnergyL2') , param_global.GEnergyL2 =1-1e-4; end
if ~isfield(param_global,'wprojection' ), param_global.wprojection = [] ; end % wprojection

% flags
if ~isfield(input_flags,'computeOperatorNorm' ), param_global.computeOperatorNorm = 1 ; end
if ~isfield(input_flags,'solveMinimization' ), param_global.solveMinimization = 1 ; end
if ~isfield(input_flags,'dr' ), param_global.dr = 0 ; end % data reduction
if ~isfield(input_flags,'homotopy' ), param_global.homotopy = 0 ; end % hootopy strategy in re-weighting


%project dir
if ~isfield(param_global,'main_dir' ), param_global.main_dir = [pwd,filesep] ; end

% default: real exp
if ~isfield(param_global,'exp_type' ), param_global.exp_type = 'realexp' ; end

% pre-processing
if ~isfield(param_global,'l2bounds_filename' ), param_global.l2bounds_filename = [] ; end
if ~isfield(param_global,'die_filename' ), param_global.die_filename = [] ; end
if ~isfield(param_global,'model_filename' ), param_global.model_filename = [] ; end
if ~isfield(param_global,'dde_filename' ), param_global.dde_filename = [] ;end


%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%% get params

% meas. op.: data blocks
% data blocks
nDataBlk= param_global.nDataBlk ;
szDataBlk= param_global.sizeDataBlk ;
l2bounds_filename = param_global.l2bounds_filename ; % available l2bounds
% meas. op.: w projection
param_wproj.CEnergyL2 = param_global.CEnergyL2 ;
param_wproj.GEnergyL2 = param_global.GEnergyL2 ;
param_wproj.do = param_global.wprojection;
% facet related
Qx = param_global.Qx;
Qy = param_global.Qy;
Qc = param_global.Qc;
window_type = param_global.window_type;
overlap_fraction = param_global.overlap_fraction;
% image dimensions & resolution
Nx = param_global.Nx;
Ny = param_global.Ny;
pixelSize  = param_global.pixelSize;
switch imageResolution
    case 'nominal'
        fprintf('\nWARNING: No pixelsize provided by user --> adopting 2x instrumental resolution.\n')
    otherwise
        fprintf('\nINFO: Pixelsize provided by user: %f asec.\n',pixelSize);
end
%* Prior
dict.nlevel = param_global.wavelet_level ; % depth of the wavelet decompositions
dict.basis  = param_global.wavelet_basis; % %! always specify Dirac basis ('self') in last position if ever used
dict.filter_length = [2*(1:(numel(dict.basis)-1))'; 0]; % length of the filters (0 corresponding to the 'self' basis)
% flags
flag_computeOperatorNorm = input_flags.computeOperatorNorm;
flag_solveMinimization = input_flags.solveMinimization;
flag_dataReduction = input_flags.dr;
flag_homotopy =input_flags.homotopy ;
exp_type = param_global.exp_type;
%% pre-processing step
param_preproc.die_filename  = param_global.die_filename;
param_preproc.l2bounds_filename = param_global.l2bounds_filename;
param_preproc.model_filename = param_global.model_filename;
param_preproc.subcube = subcube_ind;
param_preproc.done = ~isempty(param_preproc.model_filename)*( (~isempty(param_preproc.die_filename ) || ~isempty(param_preproc.dde_filename )));
%% info
format compact;
disp('General info: Experiment''s settings')
disp(['Algorithm version: ', algo_version]);
disp(['Reference image: ', image_name]);
disp(['Number of facets Qy x Qx : ', num2str(Qy), ' x ', num2str(Qx)]);
disp(['Number of spectral facets Qc : ', num2str(Qc)]);
if ~strcmp(algo_version, 'sara')
    disp(['Overlap fraction: ', strjoin(strsplit(num2str(overlap_fraction)), ', ')]);
end
%% paths
fprintf('\nSetting paths ..\n')
project_dir = param_global.main_dir;
fprintf('\nMain project dir is %s: ',project_dir)
current_dir = pwd;
fprintf('\nCurrent dir is  %s: ',current_dir)

% addpath(genpath([dir_Project,filesep,'lib']));
addpath([project_dir,filesep,'lib',filesep,'operators',filesep])
addpath([project_dir,filesep,'lib',filesep,'measurement-operator',filesep,'nufft'])
addpath([project_dir,filesep,'lib',filesep,'measurement-operator',filesep,'lib',filesep,'utils']) % added by AD
addpath([project_dir,filesep,'lib',filesep,'measurement-operator',filesep,'lib',filesep,'operators'])
% addpath ../../lib/measurement-operator/irt/nufft/
addpath([project_dir,'lib',filesep,'utils',filesep,''])
addpath([project_dir,'lib',filesep,'faceted-wavelet-transform',filesep,'src'])
addpath([project_dir,filesep,'src'])
addpath([project_dir,filesep,'src',filesep,'heuristics',filesep])
addpath([project_dir,filesep,'src',filesep,'heuristics',filesep])

if strcmp(algo_version, 'sara')
    addpath([project_dir,filesep,'src',filesep,'sara'])
elseif strcmp(algo_version, 'hs')
    addpath([project_dir,filesep,'src',filesep,'hs',filesep])
else,   addpath([project_dir,filesep,'src',filesep,'fhs',filesep])
end
data_path=fullfile(project_dir,'data',filesep);
addpath(data_path)
% setting paths to results and reference image cube
results_path = fullfile('results', strcat(image_name, '_', exp_type));
auxiliary_path = fullfile(results_path, algo_version);
mkdir(data_path)
mkdir(results_path)
mkdir(auxiliary_path)

% AD: real data utils
addpath([current_dir,filesep,'real_data'])
addpath([current_dir,filesep,'real_data',filesep,'wproj_utilities'])

%%
%! see how to load the cube prepared by Arwa
%! load the appropriate portion of the reference image cube
% if spectral faceting, just load the intersting portion of the full image cub
% AD: loading option need to be removed from the final scripts
N = Nx*Ny;
nchans  = numel(spwins2image);
ImageCubeDims = [Ny,Nx,nchans];

% load real_data;%:AD: modified path
for iSpWin = 1:nchans
    % load data and pixelsize, the rest later ..
    dataSpWinloaded=load([data_file,num2str(spwins2image(iSpWin)),'.mat'],'y', 'pos');
    fprintf('\nINFO: reading data from a saved mat file..\nINFO: uv coor. are scaled in [-pi pi].')
    y{iSpWin,1} = dataSpWinloaded.y; dataSpWinloaded.y = [];   
    pos{iSpWin,1}= dataSpWinloaded.pos; dataSpWinloaded.pos =[];
end, clear dataSpWinloaded;
% l2 bounds
if ~isempty(l2bounds_filename)
    fprintf('\nLoading estimates of the l2-bounds ...')
    for iSpWin =1:nchans
        l2SpWinloaded = load(param_preproc.l2bounds_filename(spwins2image(iSpWin)),'sigmac','l2bounds','RESIDUAL_DATA');
        for iCh=1:numel(y{iSpWin})
            for iConfig = 1 : numel(y{iSpWin}{iCh})
                try  % get residual data
                    residual_data{iSpWin}{iCh}{iConfig}  = l2SpWinloaded.RESIDUAL_DATA{iConfig,iCh};
                    l2SpWinloaded.RESIDUAL_DATA{iConfig,iCh} =[];
                    % get bounds
                    if ~flag_dataReduction, l2bounds{iSpWin}{iCh}{iConfig} = l2SpWinloaded.l2bounds.cont{iConfig,iCh};
                    else, l2bounds{iSpWin}{1}{1} =  l2SpWinloaded.l2bounds.gridded;
                    end
                catch,l2bounds{iSpWin}{iCh}{iConfig} =0;
                end
            end
        end
        try % get noise std
            if ~flag_dataReduction,	global_sigma_noise(iSpWin,1) = full(l2SpWinloaded.sigmac.cont);
            else, global_sigma_noise(iSpWin,1)=full(l2SpWinloaded.sigmac.gridded);
            end
        catch,  global_sigma_noise(iSpWin,1) =1;
        end
        clear l2SpWinloaded;
    end
else
    sigma_ball = 2;
    % assuming Chi squared dist.
    for iSpWin =  1 : nchans
        global_sigma_noise(iSpWin) =1;
        for iCh=1:numel(y{iSpWin})
            for iConfig = 1 : numel(y{iSpWin}{iCh})
                l2bounds{iSpWin}{iCh}{iConfig} =  sqrt(pos{iSpWin}{iCh}(iConfig) + sigma_ball*sqrt(pos{iSpWin}{iCh}(iConfig))) ;
            end
        end
    end
end


% get CORRECTED_DATA  if DIEs available
if ~isempty(param_preproc.die_filename )
    fprintf('\nLoading and applying DIE solutions ...')
    for iSpWin =1:nchans
        load(param_preproc.die_filename(spwins2image(iSpWin)),'DIES');
        % apply to data
        for iCh =1:numel(y{iSpWin})
            for iConfig = 1 : numel(y{iSpWin}{iCh})
                try    y{iSpWin}{iCh}{iConfig} = y{iSpWin}{iCh}{iConfig}./(DIES{iConfig,iCh});
                catch, fprintf('\nWARNING: Correcting data:  DIE solutions not read properly ')
                end
                DIES{iConfig,iCh} =[];
            end
        end
    end
end, clear DIES;


% restructure data with multiple blocks
for iSpWin =1:nchans
    y{iSpWin,1} = collapse_cell_dim(y{iSpWin,1}); %y{iSpWin,1}=vertcat(y{iSpWin,1}{:}); y{iSpWin,1} = y{iSpWin,1}(:);
    pos{iSpWin,1} = collapse_cell_dim(pos{iSpWin,1});
    l2bounds{iSpWin,1} = collapse_cell_dim(l2bounds{iSpWin,1});
end

%% image cube initialisation
switch algo_version
    case 'sara'
        xinit=zeros(Ny,Nx);
        if param_preproc.done
            try
                if isa(param_preproc.model_filename, 'function_handle')
                    xinit = fitsread(param_preproc.model_filename(spwins2image));
                else, xinit = fitsread(param_preproc.model_filename);
                end
                if  ~(size(xinit,1)== Ny && size(xinit,2)== Nx)
                    xinit=zeros(Ny,Nx);
                    fprintf('\nWARNING: init. image has a different size')
                end
            catch, fprintf('\nWARNING: init. image not found')
            end
            
        end
    otherwise
        if param_preproc.done
            % loading each channel separately 
            if isa(param_preproc.model_filename, 'function_handle')
                xinit = zeros(N, nchans);
                for  iSpWin =  1 : nchans
                    try    xinit(:,iSpWin) = reshape(fitsread(param_preproc.model_filename(spwins2image(iSpWin))),N,1);
                           fprintf('\nReading each channel initial model from a fits file.')
                    catch, fprintf('\nWARNING: init. image is of different size or not found.')
                    end
                end
            else
                % loading  a wideband cube
                try xinit = fitsread(param_preproc.model_filename);
                    xinit = reshape(xinit,[N,nchans]);
                    fprintf('\nReading the wideband initial model cube from a fits file.')
                catch, xinit = zeros(N, nchans);
                      fprintf('\nWARNING: init. image is of different size or not found.')
                end
            end
        end
end
%% Auxiliary function needed to select the appropriate workers
% (only needed for 'hs' and 'fhs' algorithms)
switch algo_version
    case 'sara'
        data_worker_id = @(k) k;
    case 'hs'
        data_worker_id = @(k) k;
    case 'fhs'
        data_worker_id = @(k) k + Qx*Qy;
    otherwise
        error('Undefined algo_version');
end
%% Get faceting parameter (spectral + spatial)
% fix faceting parameters in case these are not consistent with the
% selected algorithm
if strcmp(algo_version, 'sara')
    window_type = 'none';
    Qc = nchans;
    Qx = 1;
    Qy = 1;
elseif strcmp(algo_version, 'hs')
    window_type = 'none';
    Qx = 1;
    Qy = 1;
end
Q = Qx*Qy;

% convert fraction of overlap between consecutive facets into a number of pixels
overlap_size = get_overlap_size([Ny, Nx], [Qy, Qx], overlap_fraction);
disp(['Number of pixels in overlap: ', strjoin(strsplit(num2str(overlap_size)), ' x ')]);

% index of channels from the subcube to be handled on each data worker
if strcmp(algo_version, 'sara')
    ncores_data = numel(param_global.wavelet_basis);
end
freqRangeCores = split_range(ncores_data, nchans); %

%% Setup name of results file

if strcmp(algo_version,'sara')% AD
    results_name_function = @(nchannels) strcat(exp_type,'_',image_name,'_', ...
        algo_version, '_pixelsize', num2str(pixelSize), ...
        '_Ny',num2str(Ny), '_Nx',num2str(Nx), ...
        '_ind', num2str(subcube_ind),'_ch',num2str(spwins2image),...
        '_g', num2str(param_reg.gam), '_gb', num2str(param_reg.gam_bar),'.mat');
    
    temp_results_name = @(nchannels) strcat(exp_type,'_',image_name,'_', ...
        algo_version, '_pixelsize', num2str(pixelSize), ...
        '_Ny',num2str(Ny), '_Nx',num2str(Nx), ...
        '_ind', num2str(subcube_ind),'_ch',num2str(spwins2image),...
        '_g', num2str(param_reg.gam), '_gb', num2str(param_reg.gam_bar),'_');
    
    
else
    results_name_function = @(nchannels) strcat(exp_type,'_',image_name,'_', ...
        algo_version, '_pixelsize', num2str(pixelSize), ...
        '_Ny',num2str(Ny), '_Nx',num2str(Nx), '_L',num2str(nchannels), ...
        '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
        '_ind', num2str(subcube_ind), '_g', num2str(param_reg.gam), '_gb', num2str(param_reg.gam_bar), ...
        '_overlap', strjoin(strsplit(num2str(overlap_fraction)), '_'),'.mat');
    
    temp_results_name = @(nchannels) strcat(exp_type,'_',image_name,'_', ...
        algo_version, '_pixelsize', num2str(pixelSize), ...
        '_Ny',num2str(Ny), '_Nx',num2str(Nx), '_L',num2str(nchannels), ...
        '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
        '_ind', num2str(subcube_ind), '_g', num2str(param_reg.gam), '_gb', num2str(param_reg.gam_bar), ...
        '_overlap', strjoin(strsplit(num2str(overlap_fraction)), '_'));
    
end
warm_start = @(nchannels) strcat(temp_results_name(nchannels),'_rw', num2str(param_reg.rw), '.mat');

results_name = results_name_function(nchans);

%% Define problem configuration (rng, nufft, preconditioning, blocking,NNLS (epsilon estimation))
%init
param_nufft =[];
param_blocking =[];
param_precond =[];
param_nnls =[];
%
parameters_problem_dev_ad;
% param_blocking: config blocking already
if isempty(nDataBlk)
    param_blocking =[]; % no further blocking required
end
% nufft kernel
param_nufft.kernel = kernel;
% FoV info
param_wproj.FoVx = pixelSize * Nx *pi/180/3600;
param_wproj.FoVy = pixelSize * Ny *pi/180/3600;
param_wproj.uGridSize   = 1/(param_nufft.ox*param_wproj.FoVx);
param_wproj.vGridSize   = 1/(param_nufft.oy*param_wproj.FoVy);
disp('param_nufft')
disp(param_nufft)
disp('param_wproj')
disp(param_wproj)
disp('param_blocking')
disp(param_blocking)
disp('param_precond')
disp(param_precond)

%% setup parpool
delete(gcp('nocreate'))
if strcmp(algo_version,'sara'), cirrus =0; end
cirrus_cluster = util_set_parpool_dev(algo_version, ncores_data, Qx*Qy, cirrus);

%% Setup measurement operator
% TODO: define lambda function measurement operator
switch algo_version
    case 'sara'
        for iSpWin = 1:nchans
            dataSpWinloaded=load([data_file,num2str(spwins2image(iSpWin)),'.mat'] ,'pixelSize','u', 'v', 'w', 'nW', 'y');
            try ratio = pixelSize/dataSpWinloaded.pixelSize;%only for the CYGA data, u v w should not be normalized between [-pi pi]
            catch, ratio =1;
            end
            fprintf('\nINFO: Pixelsize ratio: %d\n',ratio);
            fprintf('\nINFO: ch %d: reading data from a saved mat file..\nINFO: uv coor. are scaled in [-pi pi].',spwins2image(iSpWin))
            nW{iSpWin,1} = dataSpWinloaded.nW;dataSpWinloaded.nW =[];
            for iCh = 1: numel(dataSpWinloaded.u)
                for iConfig = 1 : numel(dataSpWinloaded.u{iCh})
                    u{iSpWin,1}{iCh}{iConfig} = ratio*dataSpWinloaded.u{iCh}{iConfig}; dataSpWinloaded.u{iCh}{iConfig} =[];
                    v{iSpWin,1}{iCh}{iConfig} = ratio*dataSpWinloaded.v{iCh}{iConfig}; dataSpWinloaded.v{iCh}{iConfig} =[];
                end
            end
            w{iSpWin,1} =  dataSpWinloaded.w; dataSpWinloaded.w =[];
            clear dataSpWinloaded;
            % reshape cells
            if flag_dataReduction, residual_data{iSpWin,1} = collapse_cell_dim(residual_data{iSpWin,1}); 
            end
            u{iSpWin,1} = collapse_cell_dim(u{iSpWin,1}); 
            v{iSpWin,1} = collapse_cell_dim(v{iSpWin,1}); 
            w{iSpWin,1} = collapse_cell_dim(w{iSpWin,1}); 
            nW{iSpWin,1} = collapse_cell_dim(nW{iSpWin,1}); 
        end
        if flag_dataReduction
            % ! define Sigma (weight matrix involved in DR)
            % ! define G as the holographic matrix
            param_preproc.resy = residual_data;
            clear residual_data;
            param_preproc.ch = spwins2image;
            [A, At, G, W, aW,Sigma,y,noise] = util_gen_dr_measurement_operator_dev_ad(y,u,v,w, nW, ...
                param_precond, param_blocking, 1, Nx, Ny, param_nufft,param_wproj,param_preproc);
            param_preproc = rmfield(param_preproc, 'resy' ) ;
            
            try for iSpWin = 1:nchans
                    l2bounds{iSpWin}{1} = noise.l2bounds{iSpWin}{1};
                    global_sigma_noise(iSpWin,1) = noise.sigma{iSpWin};
                end
            catch, fprintf('\nCould not update l2 bounds given residual !!!! ')
            end
        else
            param_preproc.ch = spwins2image;
            [A, At, G, W, aW] = util_gen_measurement_operator_dev_ad(u,v,w, nW, ...
                param_precond, param_blocking, 1, Nx, Ny, param_nufft,param_wproj,param_preproc);
            Sigma = [];
        end
        
    otherwise % 'hs' or 'fhs'
        % create the measurement operator operator in parallel (depending on
        % the algorithm used)
        yCmpst = Composite();  resyCmpst=Composite(); Sigma = Composite();
        
        if strcmp(algo_version, 'hs')
            % send data to workers
            if flag_dataReduction
                for iLab = 1:nchans
                    local_fc = (freqRangeCores(iLab, 1):freqRangeCores(iLab, 2));
                    resyCmpst{iLab}= residual_data(local_fc);
                end, clear residual_data;
            end
            for iLab = 1:nchans
                local_fc = (freqRangeCores(iLab, 1):freqRangeCores(iLab, 2));
                yCmpst{iLab}= y(local_fc);
            end, clear y local_fc;
            
            % building measurement op.
            spmd
                local_fc = (freqRangeCores(labindex, 1):freqRangeCores(labindex, 2));
                param_preproc_cmpst = param_preproc;
                param_preproc_cmpst.ch = spwins2image(local_fc);
                for ifc  =1:numel(local_fc)
                    dataSpWinloaded=load([data_file,num2str(spwins2image(local_fc(ifc))),'.mat'] ,'pixelSize','u', 'v', 'w', 'nW', 'y');
                    try ratio =    pixelSize/dataSpWinloaded.pixelSize;%only for the CYGA data, u v w should not be normalized between [-pi pi]
                    catch, ratio =1;
                    end
                    w{ifc,1} = dataSpWinloaded.w;      dataSpWinloaded.w =[];
                    nW{ifc,1} = dataSpWinloaded.nW;    dataSpWinloaded.nW =[];
                    
                    fprintf('\nINFO: reading data from a saved mat file..\nINFO: uv coor. are scaled in [-pi pi].')
                    for iCh = 1: numel(dataSpWinloaded.u)
                        for iConfig = 1 : numel(dataSpWinloaded.u{iCh})
                            u{ifc,1}{iCh}{iConfig} = ratio*dataSpWinloaded.u{iCh}{iConfig}; dataSpWinloaded.u{iCh}{iConfig} =[];
                            v{ifc,1}{iCh}{iConfig} = ratio*dataSpWinloaded.v{iCh}{iConfig};dataSpWinloaded.v{iCh}{iConfig} =[];
                        end
                    end,  dataSpWinloaded = [];
                    
                    % reshape cells
                    u{ifc,1} = collapse_cell_dim(u{ifc,1}); 
                    v{ifc,1} = collapse_cell_dim(v{ifc,1}); 
                    w{ifc,1} = collapse_cell_dim(w{ifc,1}); 
                    nW{ifc,1} = collapse_cell_dim(nW{ifc,1}); 
                    if flag_dataReduction,resyCmpst{ifc,1} = collapse_cell_dim(resyCmpst{ifc,1}); 
                    end
                end
                
                % TODO: update util_gen_measurement_operator to enable kaiser kernels
                if flag_dataReduction
                    % ! define Sigma (weight matrix involved in DR)
                    % ! define G as the holographic matrix
                    param_preproc_cmpst.resy  = resyCmpst; resyCmpst =[];
                    [A, At, G, W, aW,Sigma,yCmpst] = util_gen_dr_measurement_operator_dev_ad(yCmpst,u, v,w,nW, ...
                        param_precond, param_blocking, numel(local_fc), Nx, Ny, param_nufft,param_wproj,param_preproc_cmpst);
                    param_preproc_cmpst = rmfield(param_preproc_cmpst,'resy');
                else
                    Sigma = [];
                    % ! ideally, simplify irt nufft interface to do so
                    [A, At, G, W, aW] = util_gen_measurement_operator_dev_ad(u, v,w,nW, ...
                        param_precond, param_blocking, numel(local_fc), Nx, Ny, param_nufft,param_wproj,param_preproc_cmpst);
                end,   param_preproc_cmpst=[]; u=[]; v=[]; w=[]; nW=[]; 
            end
        else
            
            % send data to workers
            if flag_dataReduction
                resyCmpst=Composite();
                for iLab = Q+1:Q+nchans
                    local_fc = (freqRangeCores(iLab-Q, 1):freqRangeCores(iLab-Q, 2));
                    resyCmpst{iLab}= residual_data(local_fc);
                end, clear residual_data;
            end
            
            for iLab = Q+1:Q+nchans
                local_fc = (freqRangeCores(iLab-Q, 1):freqRangeCores(iLab-Q, 2));
                yCmpst{iLab}= y(local_fc);
            end, clear local_fc y ;
            
            % building measurement op.
            spmd
                % define operator on data workers only
                if labindex > Q
                    local_fc = (freqRangeCores(labindex-Q, 1):freqRangeCores(labindex-Q, 2));
                    param_preproc_cmpst = param_preproc;
                    param_preproc_cmpst.ch = local_fc;
                    for ifc  = 1:numel(local_fc)
                        dataSpWinloaded = load([data_file,num2str(spwins2image(local_fc(ifc))),'.mat'] ,'pixelSize','u', 'v', 'w', 'nW','y');% 'time');
                        try ratio =   pixelSize/dataSpWinloaded.pixelSize;%only for the CYGA data, u v w should not be normalized between [-pi pi]
                        catch, ratio =1;
                        end
                        fprintf('\nINFO: Pixelsize %f, ratio: %f\n',pixelSize,ratio);
                        fprintf('\nINFO: reading data from a saved mat file..\nINFO: uv coor. are scaled in [-pi pi].')
                        w{ifc,1} = dataSpWinloaded.w;     dataSpWinloaded.w =[];
                        nW{ifc,1} = dataSpWinloaded.nW;   dataSpWinloaded.nW =[];
                        for iCh = 1: numel(dataSpWinloaded.u)
                            for iConfig = 1 : numel(dataSpWinloaded.u{iCh})
                                u{ifc,1}{iCh}{iConfig} = ratio*dataSpWinloaded.u{iCh}{iConfig}; dataSpWinloaded.u{iCh}{iConfig} =[];
                                v{ifc,1}{iCh}{iConfig} = ratio*dataSpWinloaded.v{iCh}{iConfig};dataSpWinloaded.v{iCh}{iConfig} =[];
                            end
                        end
                        dataSpWinloaded =[];
                        % reshape
                        u{ifc,1} = collapse_cell_dim(u{ifc,1});
                        v{ifc,1} = collapse_cell_dim(v{ifc,1});
                        w{ifc,1} = collapse_cell_dim(w{ifc,1});
                        nW{ifc,1} = collapse_cell_dim(nW{ifc,1});
                        if flag_dataReduction, resyCmpst{ifc,1} = collapse_cell_dim( resyCmpst{ifc,1});
                        end
                    end
                    
                    if flag_dataReduction
                        % ! define Sigma (weight matrix involved in DR)
                        % ! define G as the holographic matrix
                        fprintf('\nCompute H\n')
                        param_preproc_cmpst.resy   = resyCmpst; resyCmpst =[];
                        [A, At, G, W, aW,Sigma,yCmpst,noise] = util_gen_dr_measurement_operator_dev_ad(yCmpst,u, v,w,nW, ...
                            param_precond, param_blocking, numel(local_fc), Nx, Ny, param_nufft,param_wproj,param_preproc_cmpst);
                        param_preproc_cmpst = rmfield(param_preproc_cmpst,'resy');
                        
                        for  ifc  =1:numel(local_fc)
                            try l2bounds_{ifc} = noise.l2bounds{ifc};
                                global_sigma_noise_(ifc,1) = noise.sigma{ifc};
                            catch,    fprintf('\nWARNING: Could not update l2 bounds given residual data ')
                            end
                        end
                        
                    else
                        Sigma =[];
                        % ! ideally, simplify irt nufft interface to do so
                        [A, At, G, W, aW] = util_gen_measurement_operator_dev_ad(u, v,w,nW, ...
                            param_precond, param_blocking, numel(local_fc), Nx, Ny, param_nufft,param_wproj,param_preproc_cmpst);
                    end, param_preproc_cmpst=[]; u=[]; v=[]; w=[]; nW=[]; 
                end
            end
        end
end
disp('Data loaded successfully')

%% Free memory
clear residual_data u v w nW dummy;
clear dataSpWinloaded local_fc ;%% Free memory
clear param_wproj param_preproc*  param_blocking param_precond;%% Free memory

%% ids
subcube_channels = 1:nchans;
%% load l2 bounds (generate only full spectral dataset)
% only generatr data in 'hs' or 'fhs' configuration (otherwise, load the data)
% datafile = matfile(fullfile(results_path,l2bounds_filename));
switch algo_version
    case 'sara'
        % ! to be verified
        % all the variables are stored on the main process for sara
        epsilons = l2bounds(1);
    otherwise
        epsilons = Composite();
        for k = 1:ncores_data
            try
                epsilons{data_worker_id(k)} = l2bounds_{data_worker_id(k)};
                if k ==1, global_sigma_noise=global_sigma_noise_{data_worker_id(k)};
                else, global_sigma_noise=[global_sigma_noise ;global_sigma_noise_{data_worker_id(k)}];
                end
            catch,epsilons{data_worker_id(k)} = l2bounds(freqRangeCores(k, 1):freqRangeCores(k, 2));
            end
        end, clear global_sigma_noise_ l2bounds_;
end
%% Compute operator norm
if  numel(spwins2image) ==1
    subcube_ind=spwins2image; % for output filenames
else,subcube_ind =0;
end
fprintf('\nComputing measurement operator''s norm ..')
if strcmp(algo_version, 'sara')
    if flag_computeOperatorNorm
        [Anorm, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-6, 200,flag_dataReduction,Sigma);%AD: changed tolerance to 1e-6 instead of 1e-8 O
        save(fullfile(results_path, ...
            strcat('Anorm_', algo_version, ...
            '_Ny',num2str(Ny),'_Nx',num2str(Nx), ...
            '_ch', num2str(spwins2image), '.mat')), ...
            '-v7.3', 'Anorm', 'squared_operator_norm', 'rel_var', ...
            'squared_operator_norm_precond', 'rel_var_precond');
        clear rel_var
    else
        load(fullfile(results_path, ...
            strcat('Anorm_', algo_version, ...
            '_Ny',num2str(Ny),'_Nx',num2str(Nx), ...
            '_ch', num2str(spwins2image), '.mat')), ...
            'Anorm', 'squared_operator_norm_precond', 'squared_operator_norm');
    end
else
    if flag_computeOperatorNorm
        spmd
            if labindex > Qx*Qy*strcmp(algo_version, 'fhs')
                [An, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-8, 200,flag_dataReduction,Sigma);
            end
        end
        
        % save operator norm from the different subcubes into a single .mat
        % file
        opnormfile = matfile(fullfile(results_path, strcat('Anorm', ...
            '_Ny',num2str(Ny), '_Nx',num2str(Nx), ...
            '_L', num2str(nchans), '.mat')), 'Writable', true);
        
        opnormfile.squared_operator_norm = zeros(nchans, 1);
        opnormfile.rel_var = zeros(nchans, 1);
        opnormfile.squared_operator_norm_precond = zeros(nchans, 1);
        opnormfile.rel_var_precond = zeros(nchans, 1);
        
        Anorm = 0;
        for k = 1:ncores_data
            opnormfile.squared_operator_norm(freqRangeCores(k, 1):freqRangeCores(k, 2), 1) = squared_operator_norm{data_worker_id(k)};
            opnormfile.rel_var(freqRangeCores(k, 1):freqRangeCores(k, 2), 1) = rel_var{data_worker_id(k)};
            
            opnormfile.squared_operator_norm_precond(freqRangeCores(k, 1):freqRangeCores(k, 2), 1) = squared_operator_norm_precond{data_worker_id(k)};
            opnormfile.rel_var_precond(freqRangeCores(k, 1):freqRangeCores(k, 2), 1) = rel_var_precond{data_worker_id(k)};
            
            Anorm = max(Anorm, An{data_worker_id(k)});
        end
        clear An rel_var rel_var_precond squared_operator_norm_precond
        squared_operator_norm =   opnormfile.squared_operator_norm(subcube_channels, 1);
        squared_operator_norm_precond = opnormfile.squared_operator_norm_precond(subcube_channels, 1);
        
    else
        opnormfile = matfile(fullfile(results_path, strcat('Anorm', ...
            '_Ny',num2str(Ny), '_Nx',num2str(Nx), ...
            '_L', num2str(nchans), '.mat')));
        
        squared_operator_norm_precond = opnormfile.squared_operator_norm_precond(subcube_channels, 1);
        rel_var_precond = opnormfile.rel_var_precond(subcube_channels, 1);
        Anorm = max(squared_operator_norm_precond.*(1 + rel_var_precond));
        squared_operator_norm = opnormfile.squared_operator_norm(subcube_channels, 1);
    end
end

fprintf('Convergence parameter (measurement operator): %e \n', Anorm);


% TODO: to be added back if needed (new auxiliary functions needed)
% %% Generate initial epsilons by performing imaging with NNLS on each data block separately
% if generate_eps_nnls
%     % solve nnls per data block
%     for i = 1:nchans
%         eps_b{i} = cell(length(G{i}),1);
%         for j = 1 : length(G{i})
%             % printf('solving for band %i\n\n',i)
%             [~,eps_b{i}{j}] = fb_nnls_blocks(y{i}{j}, A, At, G{i}{j}, W{i}{j}, param_nnls);
%         end
%     end
%     mkdir('data/')
%     save('data/eps.mat','-v7.3', 'eps_b');
% end


%% Regularization parameters and solver

% estimate noise level (set regularization parameters to the same value)
% compute sig and sig_bar (estimate of the "noise level" in "SVD" and
% SARA space) involved in the reweighting scheme

if strcmp(algo_version, 'sara')
    
    % SARA dicionary (created out of the solver for SARA)
    dwtmode('zpd')
    [Psi1, Psit1] = op_p_sp_wlt_basis(dict.basis, dict.nlevel, Ny, Nx);
    P = numel(Psi1);
    Psi = cell(P, 1);
    Psit = cell(P, 1);
    s = zeros(P, 1); % number of wavelet coefficients for each dictionary
    
    for k = 1 : P
        f = '@(x_wave) HS_forward_sparsity(x_wave,Psi1{';
        f = sprintf('%s%i},Ny,Nx);', f,k);
        Psi{k} = eval(f);
        s(k) = size(Psit1{k}(zeros(Ny,Nx,1)),1);
        ft = ['@(x) HS_adjoint_sparsity(x,Psit1{' num2str(k) '},s(' num2str(k) '));'];
        Psit{k} = eval(ft);
    end
    
    % noise level / regularization parameter
    sig = compute_noise_level_sara(global_sigma_noise, squared_operator_norm);
    
    % apply multiplicative factor for the regularization parameter (if needed)
    mu = param_reg.gam*sig;
    fprintf('Noise level: sig = %e\n', sig);
    fprintf('Additional multiplicative regularisation factor gam = %e\n', param_reg.gam);
    fprintf('Regularization parameter mu = %e\n', mu);
    fprintf('Algo: %s, alpha = %.4e, mu = %.4e, sig = %.4e\n', algo_version, param_reg.gam, mu, sig);
end

if strcmp(algo_version, 'hs') || strcmp(algo_version, 'fhs')
    
    % noise level / regularization parameter
    [sig, sig_bar, mu_chi, sig_chi, sig_sara] = ...
        compute_noise_level(Ny, Nx, nchans, global_sigma_noise(:), ...
        algo_version, Qx, Qy, overlap_size, squared_operator_norm(:));
    
    % apply multiplicative factor for the regularization parameters (if needed)
    mu_bar = param_reg.gam_bar*sig_bar;
    mu = param_reg.gam*sig;
    fprintf('mu_chi = %.4e, sig_chi = %.4e, sig_sara = %.4e\n', mu_chi, sig_chi, sig_sara);
    fprintf('Noise levels: sig = %.4e, sig_bar = [%.4e, %.4e]\n', sig, min(sig_bar), max(sig_bar));
    fprintf('Additional multiplicative actors gam = %.4e, gam_bar = %.4e\n', param_reg.gam, param_reg.gam_bar);
    fprintf('Regularization parameters: mu = %.4e, mu_bar = %.4e\n', mu, mu_bar);
    fprintf('Algo: %s, gam = %.4e, gam_bar = %.4e, mu = %.4e, mu_bar = [%.4e, %.4e]\n', algo_version, param_reg.gam, param_reg.gam_bar, mu, min(mu_bar), max(mu_bar));
end


%% Define parameters for the solver (nReweights needed here)
parameters_solver_dev_ad


%%
% TODO: update solver interface
name_checkpoint = fullfile(auxiliary_path, temp_results_name(nchans));
name_warmstart = fullfile(auxiliary_path, warm_start(nchans));

if flag_solveMinimization
    %%
    mkdir('results/')
    
    if strcmp(algo_version, 'sara')
        disp('SARA')
        disp('-----------------------------------------')
        
        % ! in this case, ncores_data corresponds
        % ! to the number of workers for the wavelet transform (9 maximum)
        xsol = sara(y, epsilons, A, At, aW, G, W, Psi, Psit, ...
            param_solver, name_warmstart, name_checkpoint, param_reg.gam, ...
            flag_dataReduction, Sigma, xinit);
        
        
        fitswrite(xsol,fullfile(auxiliary_path, strcat('x_', image_name, '_', algo_version, ...
            '_pixelsize', num2str(pixelSize), ...  % '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
            '_ch', num2str(spwins2image), ...
            '_gam', num2str(param_reg.gam), ...
            '.fits')))
    else
        %%
        % spectral tesselation (non-overlapping)
        % ! to be updated tonight (need to be careful about the different variables needed + implicit parallelization conventions)
        cell_c_chunks = cell(ncores_data, 1); % ! to check
        for k = 1:ncores_data
            cell_c_chunks{k} = freqRangeCores(k, 1):freqRangeCores(k, 2);
        end
        
        %%
        xinit = reshape(xinit,ImageCubeDims);
        disp(name_warmstart)
        disp([name_checkpoint '_xsol.fits'])
        switch algo_version
            case 'hs'
                disp('HyperSARA')
                disp('-----------------------------------------')
                xsol = hyperSARA(yCmpst, epsilons, ...
                    A, At, aW, G, W, param_solver, ...
                    ncores_data, dict.basis, dict.nlevel, cell_c_chunks, ...
                    nchans, Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flag_dataReduction, Sigma, ...
                    xinit);%, [],X0);
            case 'fhs'
                disp('Faceted HyperSARA')
                disp('-----------------------------------------')
                xsol = facetHyperSARA(yCmpst, epsilons, ...
                    A, At, aW, G, W, param_solver, Qx, Qy, ncores_data, ...
                    dict.basis, dict.filter_length, dict.nlevel, window_type, ...
                    cell_c_chunks, nchans, overlap_size, param_reg.gam, param_reg.gam_bar, ...
                    Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flag_dataReduction, Sigma, ...
                    xinit);
            otherwise
                error('Unknown solver version.')
        end
        
        fitswrite(xsol,fullfile(auxiliary_path, strcat('x_', image_name, '_', algo_version, ...
            '_pixelsize', num2str(pixelSize), ...
            '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
            '_ind', num2str(subcube_ind), ...
            '_gam', num2str(param_reg.gam), '_gambar', num2str(param_reg.gam_bar), ...
            '_overlap', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
            '.fits')))
    end
end

%% aux functions
