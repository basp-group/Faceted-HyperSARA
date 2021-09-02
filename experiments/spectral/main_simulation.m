function main_simulation(image_name, nChannels, Qx, Qy, Qc, ...
    algo_version, window_type, ncores_data, ind, overlap_fraction, ...
    nReweights, coverage_path, gam, gam_bar, rw, exp_type, ...
    superresolution_factor, isnr, flag_generateVisibilities, ...
    flag_computeOperatorNorm, flag_solveMinimization, flagDR, ...
    flag_cirrus, flag_homotopy)
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
% ncores_data : int
%     Number of cores handlig the data fidelity terms ("data cores"). 
%     For Faceted HyperSARA, the total number of cores used is Qx*Qy + 
%     ncores_data + 1. For SARA and HyperSARA, represents the number of 
%     cores used for the parallelization.
% ind : int
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
%     Index of the reweighting step to restart from.
% exp_type : string ('spatial' or 'spectral')
%     Type of the experiment to be reproduced.
% superresolution_factor : double
%     Coverage superresolution factor.
% isnr : double
%     Input SNR used to generate the synthetic visibilities (value in dB).
% flag_generateVisibilities : bool
%     Flag specifying whether the visibilities need to be generated or 
%     loaded from an existing .mat file.
% flag_computeOperatorNorm : bool
%     Flag triggering the computation of the (preconditioned) operator 
%     norm.
% flag_solveMinimization : bool
%     Flag triggering the solver (SARA, HS or FHS).
% flagDR : bool
%     Flag to activate DR features in the definition of the measurement 
%     operator. 
% flag_cirrus : bool
%     Specify whether the solver runs on cirrus or not (for the creation of
%     the parpool).
% flag_homotopy : bool
%     Activate the homotopy strategy within the solver.
%
% Note
% ----
% DR features still need to be implemented in the main script. 
%

%% PARAMETERS FOR DEBUGGING
%
% image_name = 'W28_512'; %'cygASband_Cube_H'; %'W28_512';
% exp_type = 'local_test'; % 'spectral', 'spatial', 'test'
%
% Qx = 1; % 4
% Qy = 1; % 4
% Qc = 1;
% nReweights = 1;
% algo_version = 'fhs'; % 'fhs', 'hs', 'sara';
% window_type = 'triangular'; % 'hamming', 'pc'
% flag_generateVisibilities = 0;
% flag_computeOperatorNorm = 0;
% flag_solveMinimization = 1;
% ncores_data = 2; % number of cores assigned to the data fidelity terms (groups of channels)
% ind = 1; % index of the spectral facet to be reconstructed
% gam = 1;
% gam_bar = 1;
% coverage_path = "data/vla_7.95h_dt10s.uvw256.mat" ;%"data/msSpecs.mat"; % "data/vla_7.95h_dt10s.uvw256.mat";
%
% rw = -1;
% flag_homotopy = 0;
% overlap_fraction = 0;
% flagDR = 0;
% isnr = 50;
%
% nChannels = 20;
% flag_generateCube = 1;
% cubepath = @(nchannels) strcat(image_name, '_L', num2str(nchannels));
% cube_path = cubepath(nChannels);
% flag_generateCoverage = 0;
% flag_generateUndersampledCube = 0; % Default 15 channels cube with line emissions
% superresolution_factor = 2;
% flag_cirrus = false;
% kernel = 'minmax:tuned'; % 'kaiser' (for real data), 'minmax:tuned'
%%
format compact;

disp('MNRAS configuration')
disp(['Algorithm version: ', algo_version]);
disp(['Reference image: ', image_name]);
disp(['nchannels: ', num2str(nChannels)]);
disp(['Number of facets Qy x Qx : ', num2str(Qy), ' x ', num2str(Qx)]);
disp(['Number of spectral facets Qc : ', num2str(Qc)]);
disp(['Overlap fraction: ', strjoin(strsplit(num2str(overlap_fraction)), ', ')]);
disp(['Input SNR: ', num2str(isnr)]);
disp(['Generating visibilities: ', num2str(flag_generateVisibilities)]);

addpath ../../lib/operators/
addpath ../../lib/measurement-operator/nufft/
addpath ../../lib/measurement-operator/lib/operators/
addpath ../../lib/measurement-operator/lib/utils/
% addpath ../../lib/measurement-operator/irt/nufft/
addpath ../../lib/utils/
addpath ../../lib/faceted-wavelet-transform/src
addpath ../../data/
addpath ../../src/
addpath ../../src/heuristics/
if strcmp(algo_version, "sara")
    addpath ../../src/sara
elseif strcmp(algo_version, "hs")
    addpath ../../src/hs
else
    addpath ../../src/fhs
end

% setting paths to results and reference image cube
data_path = '../../data/';
results_path = fullfile('results', strcat(image_name, '_', exp_type));
auxiliary_path = fullfile(results_path, algo_version);
mkdir(data_path)
mkdir(results_path)
mkdir(auxiliary_path)


%%
%! see how to load the cube prepared by Arwa
%! load the appropriate portion of the reference image cube
% if spectral faceting, just load the intersting portion of the full image cube
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
        spectral_downsampling = 1;
        spatial_downsampling = 1;
    case "local_test"
        image_name = 'cygASband_Cube_256_512_100';
        spectral_downsampling = 25; % 5
        spatial_downsampling = 1;
        coverage_path = "data/vla_7.95h_dt10s.uvw256.mat";
    case "old_local_test"
        image_name = 'cubeW28';
        spectral_downsampling = 20;
        spatial_downsampling = 4;
        coverage_path = "data/vla_7.95h_dt10s.uvw256.mat";
    otherwise
        error("Unknown experiment type")
end

reference_cube_path = fullfile(data_path, strcat(image_name, '.fits'));
info        = fitsinfo(reference_cube_path);
rowend      = info.PrimaryData.Size(1);
colend      = info.PrimaryData.Size(2);
sliceend    = info.PrimaryData.Size(3);
if strcmp(algo_version, 'sara')
    x0 = fitsread(reference_cube_path, 'primary',...
              'Info', info,...
              'PixelRegion',{[1 spatial_downsampling rowend], ...
              [1 spatial_downsampling colend], ...
              ind});
else
    x0 = fitsread(reference_cube_path, 'primary',...
              'Info', info,...
              'PixelRegion',{[1 spatial_downsampling rowend], ...
              [1 spatial_downsampling colend], ...
              [1 spectral_downsampling sliceend]});
end
nChannels = floor(sliceend/spectral_downsampling);

[Ny, Nx, nchans] = size(x0);
N = Nx*Ny;
X0 = reshape(x0, [N, nchans]);
input_snr = isnr*ones(nchans, 1); % input SNR (in dB)

% frequency used to generate the reference cubes
nu0 = 2.052e9; % starting freq
dnu = 16e6;    % freq step
L = 100;       % number of channels
nu_vect = [nu0 (dnu*(1:L-1)+nu0)];
frequencies = nu_vect(1:floor(L/nChannels):end); % nu_vect(1:spectral_downsampling:end);

clear reference_cube_path info rowend colend sliceend 
clear spatial_downsampling spectral_downsampling


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

% index of the spectral channels involved in the subcube
interleaved_channels = split_range_interleaved(Qc, nChannels);
subcube_channels = interleaved_channels{ind};

% index of channels from the subcube to be handled on each data worker
rg_c = split_range(ncores_data, nchans);

% frequencies associated with the current subcube
fc = frequencies(subcube_channels);
fmax = frequencies(end);

% extract ground truth subcube
if Qc > 1 && ind > 0 && ~strcmp(algo_version, 'sara')
    x0 = x0(:,:,subcube_channels);
    nchans = size(x0,3);
    X0 = reshape(x0,Nx*Ny,nchans);
    input_snr = input_snr(subcube_channels);
end


%% Setup name of results file
data_name_function = @(nchannels) strcat('y_', ...
    exp_type,'_',image_name, '_srf=', num2str(superresolution_factor), ...
    '_Ny=',num2str(Ny),'_Nx=',num2str(Nx),'_L=', num2str(nchannels), ...
    '_snr=', num2str(isnr), ...
    '.mat');

results_name_function = @(nchannels) strcat(exp_type,'_',image_name,'_', ...
    algo_version,'_',window_type, '_srf=', num2str(superresolution_factor), ...
    '_Ny=',num2str(Ny), '_Nx=',num2str(Nx), '_L=',num2str(nchannels), ...
    '_Qy=', num2str(Qy), '_Qx=', num2str(Qx), '_Qc=', num2str(Qc), ...
    '_ind=', num2str(ind), '_g=', num2str(gam), '_gb=', num2str(gam_bar), ...
    '_overlap=', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
    '_hom=', num2str(flag_homotopy), ...
    '_snr=', num2str(isnr), ...
    '.mat');

temp_results_name = @(nchannels) strcat(exp_type,'_',image_name,'_', ...
    algo_version,'_',window_type, '_srf=', num2str(superresolution_factor), ...
    '_Ny=',num2str(Ny), '_Nx=',num2str(Nx), '_L=',num2str(nchannels), ...
    '_Qy=', num2str(Qy), '_Qx=', num2str(Qx), '_Qc=', num2str(Qc), ...
    '_ind=', num2str(ind), '_g=', num2str(gam), '_gb=', num2str(gam_bar), ...
    '_overlap=', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
    '_hom=', num2str(flag_homotopy), ...
    '_snr=', num2str(isnr));

warm_start = @(nchannels) strcat(temp_results_name(nchannels),'_rw=', num2str(rw), '.mat');

data_name = data_name_function(nChannels);
results_name = results_name_function(nChannels);


%% Define problem configuration (rng, nufft, preconditioning, blocking,
% NNLS (epsilon estimation), SARA dictionary)
parameters_problem


%% Generate/load uv-coverage
% generating u-v coverage
%! reminder uv-coverage and weighting
% https://casa.nrao.edu/Release4.1.0/doc/UserMan/UserMansu259.html
if flag_generateCoverage
    cov_type = 'vlaa';
    p = 0.5;
    dl = 1.1;
    hrs = 5;
    na = 27; % for vlaa
    M = na*(na-1)/2;
    % Fixing Mt = 0.5 N, take T = 0.5 N / M : na = 27 for vla
    T = floor(p*(Nx*Ny)/M); % should be > 1
    [u, v, ~] = generate_uv_coverage(T, hrs, dl, cov_type);
    u = u(:)*fc(1)/fmax;
    v = v(:)*fc(1)/fmax;
    fitswrite([u, v, ones(numel(u), 1)], coverage_path)
    disp(coverage_path);
else
    coverage_path
   
    % VLA configuration
    % A. 762775 -> 3
    % B. 268448 -> 2
    % C. 202957 -> 1
    % D. 47750 -> 0 

    if strcmp(exp_type, "spectral")
        load(coverage_path, 'uvw', 'obsId');
        size(uvw)
        u1 = uvw(obsId==2, 1)*fmax/speed_of_light;
        v1 = uvw(obsId==2, 2)*fmax/speed_of_light;  
        clear obsId
    else        
        % ! normalize u,v coverage w.r.t. the highest frequency (i.e., uv expressed in
        % units of the smallest wavelenght, associated with the highest frequency)
        load(coverage_path, 'uvw');
        size(uvw)
        u1 = uvw(:, 1)*fmax/speed_of_light;
        v1 = uvw(:, 2)*fmax/speed_of_light;  
%         load(coverage_path, 'uvw', 'obsId');
%         size(uvw)
%         u1 = uvw(obsId==3, 1)*fmax/speed_of_light;
%         v1 = uvw(obsId==3, 2)*fmax/speed_of_light;  
%         clear obsId
    end
    bmax = max(sqrt(u1.^2 + v1.^2));

    % cellsize = 3600*180/(superresolution_factor*2*pi*bmax); % in arcsec
    u = u1*pi/(superresolution_factor*bmax);
    v = v1*pi/(superresolution_factor*bmax);
    size(u)
    disp('Coverage loaded successfully')
    clear uvw u1 v1
end


%% setup parpool
cirrus_cluster = util_set_parpool(algo_version, ncores_data, Qx*Qy, flag_cirrus);


%% Setup measurement operator
% TODO: define lambda function measurement operator
switch algo_version
    case 'sara'
        if flagDR
            % ! define Sigma (weight matrix involved in DR)
            % ! define G as the holographic matrix
        else
            [A, At, G, W, aW] = util_gen_measurement_operator(u, v, ...
                param_precond, param_blocking, fc, fmax, Nx, Ny, ...
                param_nufft.Kx, param_nufft.Ky, param_nufft.ox, param_nufft.oy);
            Sigma = [];
        end
        
        % if ~flagDR
        %     apply_G = @(Fx, G) G * Fx;
        %     apply_Gdag = @(y, G, W) (G') * y(W);
        % else
        %     % ! in this case, the variable T (weights, ...) needs to be defined
        %     apply_G = @(Fx, G) T.* (G * Fx);
        %     apply_Gdag = @(y, G) G' * (T.*y);
        % end

    otherwise % 'hs' or 'fhs'

        % create the measurement operator operator in parallel (depending on
        % the algorithm used)
        if strcmp(algo_version, 'hs')
            spmd
                local_fc = fc(rg_c(labindex, 1):rg_c(labindex, 2));
                % TODO: update util_gen_measurement_operator to enable kaiser kernels
                if flagDR
                    % ! define Sigma (weight matrix involved in DR)
                    % ! define G as the holographic matrix
                else
                    % ! ideally, simplify irt nufft interface to do so
                    [A, At, G, W, aW] = util_gen_measurement_operator(u, v, ...
                    param_precond, param_blocking, local_fc, fmax, Nx, Ny, param_nufft.Kx, param_nufft.Ky, param_nufft.ox, param_nufft.oy, kernel);
                    Sigma = [];
                end
            end
        else
            spmd
                % define operator on data workers only
                if labindex > Q
                    local_fc = fc(rg_c(labindex-Q, 1):rg_c(labindex-Q, 2));
                    if flagDR
                        % ! define Sigma (weight matrix involved in DR)
                        % ! define G as the holographic matrix
                    else
                        % ! ideally, simplify irt nufft interface to do so
                        [A, At, G, W, aW] = util_gen_measurement_operator(u, v, ...
                        param_precond, param_blocking, local_fc, fmax, Nx, Ny, param_nufft.Kx, param_nufft.Ky, param_nufft.ox, param_nufft.oy, kernel);
                        Sigma = [];
                    end
                end
            end
        end
        clear local_fc
end


%% Free memory
clear param_blocking param_precond;

%% Generate/load visibilities (generate only full spectral dataset)
% only generatr data in 'hs' or 'fhs' configuration (otherwise, load the data)
if flag_generateVisibilities

    param_l2_ball.type = 'sigma';
    param_l2_ball.sigma_ball = 2;

    % TODO: modify data generation to allow reproducible parallel rng streams
    % ! not implemented/activated yet
    % https://fr.mathworks.com/help/matlab/math/creating-and-controlling-a-random-number-stream.html?searchHighlight=random%20number%20streams&s_tid=srchtitle#brvku_2
    % may not be strictly equivalent to the data initially obtained for the
    % mnras paper

    % rng_stream = RandStream.create('threefry4x64_20', ...
    %     'Seed', seed, 'NumStreams', ncores_data);
    % offset_worker = Q*strcmp(algo_version, 'fhs');

    spmd
        if labindex > Q*strcmp(algo_version, 'fhs')
            [y0, y, Ml, ~, sigma_noise, ~] = util_gen_measurements_snr( ...
                x0(:,:,rg_c(labindex, 1):rg_c(labindex, 2)), G, W, A, ...
                input_snr(rg_c(labindex, 1):rg_c(labindex, 2)));
                % rng_stream(labindex-offset_worker)
            [~, epsilons] = util_gen_data_fidelity_bounds2(y, Ml, .../
                param_l2_ball, sigma_noise);
        end
    end
    
    % save parameters (matfile solution)
    datafile = matfile(fullfile(results_path,data_name), ...
        'Writable', true);
    datafile.y0 = cell(nchans, 1);
    datafile.y = cell(nchans, 1);
    datafile.epsilons = cell(nchans, 1);
    datafile.sigma_noise = zeros(nchans, 1);

    for k = 1:ncores_data
        % convert from local to global channel index (i.e., index into the full spectral dimension using "subcube_channels(rg_c(k, 1))"
        datafile.y0(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)),1) = y0{data_worker_id(k)};
        datafile.y(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)),1) = y{data_worker_id(k)};
        datafile.epsilons(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)),1) = epsilons{data_worker_id(k)};
        datafile.sigma_noise(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1) = sigma_noise{data_worker_id(k)};
    end
    global_sigma_noise = datafile.sigma_noise;
    clear param_l2_ball m Ml epsilons datafile
else
    datafile = matfile(fullfile(results_path,data_name));
    
    switch algo_version
        case 'sara'
            % ! to be verified
            % all the variables are stored on the main process for sara
            y = datafile.y(subcube_channels, 1); % subcube_channels contains a single index for SARA
            epsilons = datafile.epsilons(subcube_channels, 1);
            global_sigma_noise = datafile.sigma_noise(subcube_channels, 1);
        otherwise
            y = Composite();
            epsilons = Composite();
            sigma_noise = Composite();

            for k = 1:ncores_data
                y{data_worker_id(k)} = datafile.y(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1);
                epsilons{data_worker_id(k)} = datafile.epsilons(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1);
                sigma_noise{data_worker_id(k)} = datafile.sigma_noise(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1);
            end
            global_sigma_noise = datafile.sigma_noise;
    end
    disp('Data loaded successfully')
end


%% Compute operator norm
if strcmp(algo_version, 'sara')
    if flag_computeOperatorNorm
        [Anorm, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-8, 200);
        
        save(fullfile(results_path, ...
            strcat('Anorm_', algo_version, ...
            '_Ny=',num2str(Ny),'_Nx=',num2str(Nx), ...
            '_L=',num2str(nChannels), ...
            '_Qc=',num2str(Qc),'_ind=',num2str(ind), ...
            '_ch=', num2str(ind), '.mat')), ...
            '-v7.3', 'Anorm', 'squared_operator_norm', 'rel_var', ...
            'squared_operator_norm_precond', 'rel_var_precond');
        clear rel_var
    else
        load(fullfile(results_path, ...
            strcat('Anorm_', algo_version, ...
            '_Ny=',num2str(Ny),'_Nx=',num2str(Nx), ...
            '_L=',num2str(nChannels), ...
            '_Qc=',num2str(Qc),'_ind=',num2str(ind), ...
            '_ch=', num2str(ind), '.mat')), ...
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
            '_Ny=',num2str(Ny), '_Nx=',num2str(Nx), ...
            '_L=', num2str(nChannels), '.mat')), 'Writable', true);

        opnormfile.squared_operator_norm = zeros(nchans, 1);
        opnormfile.rel_var = zeros(nchans, 1);
        opnormfile.squared_operator_norm_precond = zeros(nchans, 1);
        opnormfile.rel_var_precond = zeros(nchans, 1);

        Anorm = 0;
        for k = 1:ncores_data
            opnormfile.squared_operator_norm(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1) = squared_operator_norm{data_worker_id(k)};
            opnormfile.rel_var(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1) = rel_var{data_worker_id(k)};

            opnormfile.squared_operator_norm_precond(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1) = squared_operator_norm_precond{data_worker_id(k)};
            opnormfile.rel_var_precond(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1) = rel_var_precond{data_worker_id(k)};

            Anorm = max(Anorm, An{data_worker_id(k)});
        end
        clear An rel_var rel_var_precond squared_operator_norm_precond

    else
        opnormfile = matfile(fullfile(results_path, strcat('Anorm', ...
            '_Ny=',num2str(Ny), '_Nx=',num2str(Nx), ...
            '_L=', num2str(nChannels), '.mat')));

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
    [Psi1, Psit1] = op_p_sp_wlt_basis_fhs(wlt_basis, nlevel, Ny, Nx);
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
    mu = gam*sig;
    fprintf('Noise level: sig = %e\n', sig);
    fprintf('Additional multiplicative regularisation factor gam = %e\n', gam);
    fprintf('Regularization parameter mu = %e\n', mu);
    fprintf('Algo: %s, alpha = %.4e, mu = %.4e, sig = %.4e\n', algo_version, gam, mu, sig);
end

if strcmp(algo_version, 'hs') || strcmp(algo_version, 'fhs')

    % noise level / regularization parameter
    [sig, sig_bar, mu_chi, sig_chi, sig_sara] = ...
        compute_noise_level(Ny, Nx, nchans, global_sigma_noise, ...
        algo_version, Qx, Qy, overlap_size, squared_operator_norm);

    % apply multiplicative factor for the regularization parameters (if needed)
    mu_bar = gam_bar*sig_bar;
    mu = gam*sig;
    fprintf('mu_chi = %.4e, sig_chi = %.4e, sig_sara = %.4e\n', mu_chi, sig_chi, sig_sara);
    fprintf('Noise levels: sig = %.4e, sig_bar = [%.4e, %.4e]\n', sig, min(sig_bar), max(sig_bar));
    fprintf('Additional multiplicative actors gam = %.4e, gam_bar = %.4e\n', gam, gam_bar);
    fprintf('Regularization parameters: mu = %.4e, mu_bar = %.4e\n', mu, mu_bar);
    fprintf('Algo: %s, gam = %.4e, gam_bar = %.4e, mu = %.4e, mu_bar = [%.4e, %.4e]\n', algo_version, gam, gam_bar, mu, min(mu_bar), max(mu_bar));
end


%% Define parameters for the solver (nReweights needed here)
parameters_solver 


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
                    param_solver, name_warmstart, name_checkpoint, gam, ...
                    flagDR, Sigma, [], x0);

        mkdir('results/')

        fitswrite(xsol,fullfile(auxiliary_path, strcat('x_', image_name, '_', algo_version, ...
        '_srf=', num2str(superresolution_factor), ...
        '_', window_type, ...
        '_Qy=', num2str(Qy), '_Qx=', num2str(Qx), '_Qc=', num2str(Qc), ...
        '_ind=', num2str(ind), ...
        '_gam=', num2str(gam), ...
        '_homotopy=', num2str(flag_homotopy), ...
        '_snr=', num2str(isnr), ...
        '.fits')))
    else
        %%

        % spectral tesselation (non-overlapping)
        % ! to be updated tonight (need to be careful about the different variables needed + implicit parallelization conventions)
        cell_c_chunks = cell(ncores_data, 1); % ! to check
        for k = 1:ncores_data
            cell_c_chunks{k} = rg_c(k, 1):rg_c(k, 2);
        end
        
        %%
        switch algo_version 
            case 'hs'
                disp('HyperSARA')
                disp('-----------------------------------------')
                xsol = hyperSARA(y, epsilons, ...
                    A, At, aW, G, W, param_solver, ...
                    ncores_data, wlt_basis, nlevel, cell_c_chunks, ...
                    nchans, Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flagDR, Sigma, ...
                    [], X0);
            case 'fhs'
                disp('Faceted HyperSARA')
                disp('-----------------------------------------')
                xsol = facetHyperSARA(y, epsilons, ...
                    A, At, aW, G, W, param_solver, Qx, Qy, ncores_data, ...
                    wlt_basis, filter_length, nlevel, window_type, ...
                    cell_c_chunks, nchans, overlap_size, gam, gam_bar, ...
                    Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, flagDR, Sigma, ...
                    [], X0);
            otherwise
                error('Unknown solver version.')
        end

        mkdir('results/')
        fitswrite(xsol,fullfile(auxiliary_path, strcat('x_', image_name, '_', algo_version, ...
            '_', window_type, ...
            '_srf=', num2str(superresolution_factor), ...
            '_Qy=', num2str(Qy), '_Qx=', num2str(Qx), '_Qc=', num2str(Qc), ...
            '_ind=', num2str(ind), ...
            '_gam=', num2str(gam), '_gambar=', num2str(gam_bar), ...
            '_overlap=', strjoin(strsplit(num2str(overlap_fraction)), '_'), ...
            '_homotopy=', num2str(flag_homotopy), ...
            '_snr=', num2str(isnr), ...
            '.fits')))
    end       
end
