function main_real_data_dev_ad(image_name, data_file,subcube_ind,channels2image, ...
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
%     Number of spatial facets along axis 1 (y).
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

% flags

if ~isfield(input_flags,'computeOperatorNorm' ), input_flags.computeOperatorNorm = 1 ; end
if ~isfield(input_flags,'solveMinimization' ), input_flags.solveMinimization = 1 ; end
if ~isfield(input_flags,'dr' ), input_flags.dr = 0 ; end
if ~isfield(input_flags,'homotopy' ), input_flags.homotopy = 0 ; end

%project dir
if ~isfield(param_global,'main_dir' ), param_global.main_dir = [pwd,filesep] ; end

% wprojection
if ~isfield(param_global,'wprojection' ), param_global.wprojection = [] ; end
% default: real exp
if ~isfield(param_global,'exp_type' ), param_global.exp_type = 'realexp' ; end


if ~isfield(param_global,'l2bounds_filename' ), param_global.l2bounds_filename = [] ; end
if ~isfield(param_global,'G_filename' ), param_global.G_filename = [] ; end
if ~isfield(param_global,'model_filename' ), param_global.model_filename = [] ; end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% get flags
% meas. op.: data blocks
% data blocks
nDataBlk= param_global.nDataBlk ;
szDataBlk= param_global.sizeDataBlk ;
l2bounds_filename = param_global.l2bounds_filename ; % available l2bounds
% meas. op.: w projection
param_wproj.CEnergyL2 = param_global.CEnergyL2 ;
param_wproj.GEnergyL2 = param_global.GEnergyL2 ;
param_wproj.do = param_global.wprojection;
%% get params
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
flag_dr = input_flags.dr;
flag_homotopy =input_flags.homotopy ;
exp_type = param_global.exp_type;
%% pre-processing step
param_preproc.G_filename  = param_global.G_filename;
param_preproc.l2bounds_filename = param_global.l2bounds_filename;
param_preproc.model_filename = param_global.model_filename;
param_preproc.subcube = subcube_ind;
param_preproc.done = ~isempty(param_preproc.l2bounds_filename ) * ~isempty(param_preproc.model_filename )* ~isempty(param_preproc.G_filename );
%% info + paths
format compact;
project_dir = param_global.main_dir;
fprintf('\nMain project dir is %s: ',project_dir)
current_dir = pwd;
fprintf('\nCurrent dir is  %s: ',current_dir)

disp('MNRAS configuration')
disp(['Algorithm version: ', algo_version]);
disp(['Reference image: ', image_name]);
% disp(['nchannels: ', num2str(nChannels)]);
disp(['Number of facets Qy x Qx : ', num2str(Qy), ' x ', num2str(Qx)]);
disp(['Number of spectral facets Qc : ', num2str(Qc)]);
disp(['Overlap fraction: ', strjoin(strsplit(num2str(overlap_fraction)), ', ')]);

% addpath(genpath([dir_Project,filesep,'lib']));
addpath([project_dir,filesep,'lib',filesep,'operators',filesep])
addpath([project_dir,filesep,'lib',filesep,'measurement-operator',filesep,'nufft'])
addpath([project_dir,filesep,'lib',filesep,'measurement-operator',filesep,'lib',filesep,'operators'])
% addpath ../../lib/measurement-operator/irt/nufft/
addpath([project_dir,'lib',filesep,'utils',filesep,''])
addpath([project_dir,'lib',filesep,'faceted-wavelet-transform',filesep,'src'])
addpath([project_dir,filesep,'src'])
addpath([project_dir,filesep,'src',filesep,'heuristics',filesep])
addpath([project_dir,filesep,'src',filesep,'heuristics',filesep])

if strcmp(algo_version, "sara")
    addpath([project_dir,filesep,'src',filesep,'sara'])
elseif strcmp(algo_version, "hs")
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

N = Nx*Ny;
% load real_data;%:AD: modified path
load([data_file,num2str(subcube_ind),'.mat'],'y', 'u', 'v', 'w', 'nW', 'time', 'pos','pixelSize');
fprintf('\nINFO: reading data from a saved mat file..\nINFO: uv coor. are scaled in [-pi pi].')
fprintf('\nINFO: pixelsize: %f asec\n',pixelSize)

y = y(channels2image);
u = u(channels2image);
v = v(channels2image);
w = w(channels2image);
nW = nW(channels2image);
time = time(channels2image);
pos = pos(channels2image);
nchans  = numel(channels2image);
if ~isempty(l2bounds_filename)
    fprintf('\nLoading estimates of the l2-bounds ...')
    for iCh =1:nchans
        loaded = load(param_preproc.l2bounds_filename(subcube_ind,channels2image(iCh)),'l2bounds');
        l2bounds{iCh} = loaded.l2bounds;
    end
else
    sigma_ball = 2;
    % assuming Chi squared dist.
    for iCh =  1 : nchans
        %         l2bounds{iCh} = cell(numel(pos{iCh}), 1);
        for iConfig = 1 : numel(pos{iCh})
            l2bounds{iCh}(iConfig,1) =  sqrt(pos{iCh}(iConfig) + sigma_ball*sqrt(pos{iCh}(iConfig))) ;
        end
    end
end
%
nChannels =nchans;
%% image cube initialisation
switch algo_version
    case 'sara'
        if ~param_preproc.done
            xinit=zeros(Ny,Nx);
        else
            xinit =fitsread(param_preproc.model_filename(subcube_ind,channels2image));
        end
        x0 = sparse(Ny,Nx);
    otherwise
        xinit = zeros(N, nchans);
        if param_preproc.done
            for  iCh =  1 : nchans
                xinit(:,iCh) = reshape(fitsread(param_preproc.model_filename(subcube_ind,channels2image(iCh))),N,1);
            end
        end
end

X0 = sparse(N, nchans); %init cube

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
    Qc = nChannels;
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
ncores_data = min(ncores_data,nchans);
freqRangeCores = split_range(ncores_data, nchans); %

%% Setup name of results file
results_name_function = @(nchannels) strcat(exp_type,'_',image_name,'_', ...
    algo_version,'_',window_type, '_pixelsize', num2str(pixelSize), ...
    '_Ny',num2str(Ny), '_Nx',num2str(Nx), '_L',num2str(nchannels), ...
    '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
    '_ind', num2str(subcube_ind), '_g', num2str(param_reg.gam), '_gb', num2str(param_reg.gam_bar), ...
    '_overlap', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
    '_hom', num2str(flag_homotopy),'.mat');

temp_results_name = @(nchannels) strcat(exp_type,'_',image_name,'_', ...
    algo_version,'_',window_type, '_pixelsize', num2str(pixelSize), ...
    '_Ny',num2str(Ny), '_Nx',num2str(Nx), '_L',num2str(nchannels), ...
    '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
    '_ind', num2str(subcube_ind), '_g', num2str(param_reg.gam), '_gb', num2str(param_reg.gam_bar), ...
    '_overlap', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
    '_hom', num2str(flag_homotopy));

warm_start = @(nchannels) strcat(temp_results_name(nchannels),'_rw', num2str(param_reg.rw), '.mat');

results_name = results_name_function(nChannels);

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
cirrus_cluster = util_set_parpool(algo_version, ncores_data, Qx*Qy, cirrus);

%% Setup measurement operator
% TODO: define lambda function measurement operator
switch algo_version
    case 'sara'
        if flag_dr
            % ! define Sigma (weight matrix involved in DR)
            % ! define G as the holographic matrix
        else
            param_preproc.ch = channels2image;
            [A, At, G, W, aW] = util_gen_measurement_operator_dev_ad(u,v,w, nW, ...
                param_precond, param_blocking, 1, Nx, Ny, param_nufft,param_wproj,param_preproc);
            Sigma = [];
        end
        
    otherwise % 'hs' or 'fhs'
        % create the measurement operator operator in parallel (depending on
        % the algorithm used)
        
        if strcmp(algo_version, 'hs')
            spmd
                local_fc = (freqRangeCores(labindex, 1):freqRangeCores(labindex, 2));
                param_preproc.ch = channels2image(local_fc);
                
                % TODO: update util_gen_measurement_operator to enable kaiser kernels
                if flag_dr
                    % ! define Sigma (weight matrix involved in DR)
                    % ! define G as the holographic matrix
                else
                    % ! ideally, simplify irt nufft interface to do so
                    [A, At, G, W, aW] = util_gen_measurement_operator_dev_ad(u(local_fc), v(local_fc),w(local_fc),nW(local_fc), ...
                        param_precond, param_blocking, numel(local_fc), Nx, Ny, param_nufft,param_wproj,param_preproc);
                    Sigma = [];
                end
            end
        else
            spmd
                % define operator on data workers only
                if labindex > Q
                    
                    local_fc = (freqRangeCores(labindex-Q, 1):freqRangeCores(labindex-Q, 2));
                    param_preproc.ch = local_fc;
                    if flag_dr
                        % ! define Sigma (weight matrix involved in DR)
                        % ! define G as the holographic matrix
                    else
                        % ! ideally, simplify irt nufft interface to do so
                        [A, At, G, W, aW] = util_gen_measurement_operator_dev_ad(u(local_fc), v(local_fc),w(local_fc),nW(local_fc), ...
                            param_precond, param_blocking, numel(local_fc), Nx, Ny, param_nufft,param_wproj,param_preproc);
                        Sigma = [];
                    end
                end
            end
        end, clear local_fc
        
end, clear param_blocking param_precond;%% Free memory

% ids
subcube_channels = 1:nchans;
%% load visibilities (generate only full spectral dataset)
% only generatr data in 'hs' or 'fhs' configuration (otherwise, load the data)
% datafile = matfile(fullfile(results_path,l2bounds_filename));

switch algo_version
    case 'sara'
        % ! to be verified
        % all the variables are stored on the main process for sara
        epsilons = l2bounds(1);
        global_sigma_noise = 1;% data are whitened already%datafile.sigma_noise(subcube_channels, 1);
    otherwise
        yCmpst = Composite();
        epsilons = Composite();
        
        for k = 1:ncores_data
            yCmpst{data_worker_id(k)} = y(freqRangeCores(k, 1):freqRangeCores(k, 2));
            epsilons{data_worker_id(k)} = l2bounds(freqRangeCores(k, 1):freqRangeCores(k, 2));
        end
        global_sigma_noise = 1;%data are whitened already %datafile.sigma_noise;
end
disp('Data loaded successfully')



%% Compute operator norm
if strcmp(algo_version, 'sara')
    if flag_computeOperatorNorm
        [Anorm, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-6, 200);%AD: changed tolerance to 1e-6 instead of 1e-8 Oo
        save(fullfile(results_path, ...
            strcat('Anorm_', algo_version, ...
            '_Ny',num2str(Ny),'_Nx',num2str(Nx), ...
            '_L',num2str(nChannels), ...
            '_Qc',num2str(Qc),'_ind',num2str(subcube_ind), ...
            '_ch', num2str(subcube_ind), '.mat')), ...
            '-v7.3', 'Anorm', 'squared_operator_norm', 'rel_var', ...
            'squared_operator_norm_precond', 'rel_var_precond');
        clear rel_var
    else
        load(fullfile(results_path, ...
            strcat('Anorm_', algo_version, ...
            '_Ny',num2str(Ny),'_Nx',num2str(Nx), ...
            '_L',num2str(nChannels), ...
            '_Qc',num2str(Qc),'_ind',num2str(subcube_ind), ...
            '_ch', num2str(subcube_ind), '.mat')), ...
            'Anorm', 'squared_operator_norm_precond', 'squared_operator_norm');
    end
else
    if flag_computeOperatorNorm
        spmd
            if labindex > Qx*Qy*strcmp(algo_version, 'fhs')
                [An, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-8, 200);
            end
        end
        
        % save operator norm from the different subcubes into a single .mat
        % file
        opnormfile = matfile(fullfile(results_path, strcat('Anorm', ...
            '_Ny',num2str(Ny), '_Nx',num2str(Nx), ...
            '_L', num2str(nChannels), '.mat')), 'Writable', true);
        
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
        
    else
        opnormfile = matfile(fullfile(results_path, strcat('Anorm', ...
            '_Ny',num2str(Ny), '_Nx',num2str(Nx), ...
            '_L', num2str(nChannels), '.mat')));
        
        squared_operator_norm_precond = opnormfile.squared_operator_norm_precond(subcube_channels, 1);
        rel_var_precond = opnormfile.rel_var_precond(subcube_channels, 1);
        Anorm = max(squared_operator_norm_precond.*(1 + rel_var_precond));
        squared_operator_norm = opnormfile.squared_operator_norm(subcube_channels, 1);
        
        % squared_operator_norm = Composite();
        % for k = 1:ncores_data
        %     squared_operator_norm{data_worker_id(k)} = opnormfile.squared_operator_norm(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1);
        % end
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
    [Psi1, Psit1] = op_p_sp_wlt_basis_fhs(dict.basis, dict.nlevel, Ny, Nx);
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
        compute_noise_level(Ny, Nx, nchans, global_sigma_noise, ...
        algo_version, Qx, Qy, overlap_size, squared_operator_norm);
    
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
name_checkpoint = fullfile(auxiliary_path, temp_results_name(nChannels));
name_warmstart = fullfile(auxiliary_path, warm_start(nChannels));

if flag_solveMinimization
    %%
    if strcmp(algo_version, 'sara')
        disp('SARA')
        disp('-----------------------------------------')
        
        % ! in this case, ncores_data corresponds
        % ! to the number of workers for the wavelet transform (9 maximum)
        xsol = sara(y, epsilons, A, At, aW, G, W, Psi, Psit, ...
            param_solver, name_warmstart, name_checkpoint, param_reg.gam, ...
            flag_dr, Sigma, xinit);%,[], x0);
        
        mkdir('results/')
        
        fitswrite(xsol,fullfile(auxiliary_path, strcat('x_', image_name, '_', algo_version, ...
            '_pixelsize', num2str(pixelSize), ...
            '_', window_type, ...
            '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
            '_ind', num2str(subcube_ind), ...
            '_gam', num2str(param_reg.gam), ...
            '_homotopy', num2str(flag_homotopy), ...
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
        switch algo_version
            case 'hs'
                disp('HyperSARA')
                disp('-----------------------------------------')
                xsol = hyperSARA(yCmpst, epsilons, ...
                    A, At, aW, G, W, param_solver, ...
                    ncores_data, dict.basis, dict.nlevel, cell_c_chunks, ...
                    nchans, Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flag_dr, Sigma, ...
                    xinit);%, [],X0);
            case 'fhs'
                disp('Faceted HyperSARA')
                disp('-----------------------------------------')
                xsol = facetHyperSARA(yCmpst, epsilons, ...
                    A, At, aW, G, W, param_solver, Qx, Qy, ncores_data, ...
                    dict.basis, dict.filter_length, dict.nlevel, window_type, ...
                    cell_c_chunks, nchans, overlap_size, param_reg.gam, param_reg.gam_bar, ...
                    Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flag_dr, Sigma, ...
                    xinit);%, [],X0);
            otherwise
                error('Unknown solver version.')
        end
        
        mkdir('results/')
        fitswrite(xsol,fullfile(auxiliary_path, strcat('x_', image_name, '_', algo_version, ...
            '_', window_type, ...
            '_pixelsize', num2str(pixelSize), ...
            '_Qy', num2str(Qy), '_Qx', num2str(Qx), '_Qc', num2str(Qc), ...
            '_ind', num2str(subcube_ind), ...
            '_gam', num2str(param_reg.gam), '_gambar', num2str(param_reg.gam_bar), ...
            '_overlap', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
            '_homotopy', num2str(flag_homotopy), ...
            '.fits')))
    end
end
