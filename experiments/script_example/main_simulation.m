function main_simulation_draft(json_filename, ind)
% Main script to run the faceted HyperSARA approach on synthetic data.
%
% This script generates synthetic data and runs the SARA, HyperSARA or
% faceted HyperSARA approach to reconstruct an :math:`N \times L` wideband
% image cube.
%
% Parameters
% ----------
% param_general.image_name : string
%     Name of the reference synthetic image (from the data/ folder).
% n_channels : int
%     Number of spectral channels considered.
% param_model.Qx : int
%     Number of spatial facets along axis 2 (x).
% param_model.Qy : int
%     Number of spatial facets along axis 1 (y).
% param_model.Qc : int
%     Number of spectral facets.
% param_model.algo_version : string ('sara', 'hs' or 'fhs')
%     Selected solver.
% param_model.window_type : string ('triangular', 'hamming' or 'pc' (piecewise-constant))
%     Type of apodization window considered for the faceted nuclear norm
%     prior (FHS solver).
% param_model.ncores_data : int
%     Number of cores handlig the data fidelity terms ("data cores").
%     For Faceted HyperSARA, the total number of cores used is param_model.Qx*param_model.Qy +
%     param_model.ncores_data + 1. For SARA and HyperSARA, represents the number of
%     cores used for the parallelization.
% ind : int
%     Index of the spectral facet to be reconstructed (set to -1 to
%     deactivate spectral faceting).
% param_model.overlap_fraction : array (1d)
%     Fraction of the total size of a facet overlapping with a neighbour
%     facet.
% param_reweighting.max_iter : int
%     Maximum number of reweighting steps.
% coverage_path : string
%     Path and name of the uv-coverage .fits file (w/o file extension).
% param_model.gamma : double
%     Additional multiplicative factor affecting the joint-sparsity
%     regularization term.
% param_model.gamma_bar : double
%     Additional multiplicative factor affecting the low-rankness
%     regularization term.
% param_model.warmstart_iteration : int
%     Index of the reweighting step to restart from.
% param_synth.exp_type : string ('spatial' or 'spectral')
%     Type of the experiment to be reproduced.
% param_synth.superresolution_factor : double
%     Coverage superresolution factor.
% param_general.flag_compute_operator_norm : bool
%     Flag triggering the computation of the (preconditioned) operator
%     norm.
% param_general.flag_solve_minimization : bool
%     Flag triggering the solver (SARA, HS or FHS).
% param_model.flagDR : bool
%     Flag to activate DR features in the definition of the measurement
%     operator.
% param_general.flag_cirrus : bool
%     Specify whether the solver runs on cirrus or not (for the creation of
%     the parpool).
%
% Note
% ----
% DR features still need to be implemented in the main script.
%

%% PARAMETERS FOR DEBUGGING
%
% param_general.image_name = 'W28_512'; %'cygASband_Cube_H'; %'W28_512';
% param_synth.exp_type = 'local_test'; % 'spectral', 'spatial', 'test'

% param_model.Qx = 1; % 4
% param_model.Qy = 1; % 4
% param_model.Qc = 1;
% param_reweighting.max_iter = 1;
% param_model.algo_version = 'fhs'; % 'fhs', 'hs', 'sara';
% param_model.window_type = 'triangular'; % 'hamming', 'pc'
% param_general.flag_compute_operator_norm = 0;
% param_general.flag_solve_minimization = 1;
% param_synth.flag_generateCoverage = 0;
% param_model.ncores_data = 2; % number of cores assigned to the data fidelity terms (groups of channels)
% ind = 1; % index of the spectral facet to be reconstructed
% param_model.gamma = 1;
% param_model.gamma_bar = 1;
% coverage_path = "data/vla_7.95h_dt10s.uvw256.mat" ;%"data/msSpecs.mat"; % "data/vla_7.95h_dt10s.uvw256.mat";

% param_model.warmstart_iteration = -1;
% param_model.overlap_fraction = 0.25;
% param_model.flagDR = 0;

% n_channels = 20;
% flag_generateCube = 1;
% cubepath = @(nchannels) strcat(param_general.image_name, '_L', num2str(nchannels));
% cube_path = cubepath(n_channels);
% param_synth.superresolution_factor = 2; % to be loaded from the data
% param_general.flag_cirrus = false;
% param_nufft.kernel = 'minmax:tuned'; % 'kaiser' (for real data), 'minmax:tuned'

%%
% TODO: update util_gen_measurement_operator to enable kaiser kernels
% json_filename = "setup_matlab.json";
% ind = 1;
format compact;

% Read .json configuration file
[speed_of_light, param_general, param_model, param_solver, param_nufft, ...
 param_blocking, param_precond, param_synth, param_nnls] = ...
    read_json_configuration(json_filename, ind);

%%
% TODO: update util_gen_measurement_operator to enable kaiser kernels

addpath ../../lib/operators/;
addpath ../../lib/measurement-operator/nufft/;
addpath ../../lib/measurement-operator/lib/operators/;
addpath ../../lib/measurement-operator/lib/utils/;
addpath ../../lib/utils/;
addpath ../../lib/faceted-wavelet-transform/src;
addpath ../../data/;
addpath ../../src/;
addpath ../../src/heuristics/;
if strcmp(param_model.algo_version, "sara")
    addpath ../../src/sara;
elseif strcmp(param_model.algo_version, "hs")
    addpath ../../src/hs;
else
    addpath ../../src/fhs;
end

% setting paths to results and reference image cube
data_path = '../../data/';
results_path = fullfile('results', strcat(param_general.image_name, '_', param_synth.exp_type));
auxiliary_path = fullfile(results_path, param_model.algo_version);
mkdir(data_path);
mkdir(results_path);
mkdir(auxiliary_path);

%%
% ! see how to load the cube prepared by Arwa
% ! load the appropriate portion of the reference image cube
% if spectral faceting, just load the intersting portion of the full image cube
switch param_synth.exp_type
    case "spatial"
        param_general.image_name = 'cygASband_Cube_1024_2048_20';
        spectral_downsampling = 1;
        spatial_downsampling = 1;
    case "spectral"
        param_general.image_name = 'cygASband_Cube_256_512_100';
        spectral_downsampling = 1;
        spatial_downsampling = 1;
    case "test"
        param_general.image_name = 'cygASband_Cube_512_1024_20';
        spectral_downsampling = 1;
        spatial_downsampling = 1;
    case "local_test"
        param_general.image_name = 'cygASband_Cube_256_512_100';
        spectral_downsampling = 25;
        spatial_downsampling = 1;
        coverage_path = "data/vla_7.95h_dt10s.uvw256.mat";
        param_model.ncores_data = 2;
        param_model.Qx = 1;
        param_model.Qy = 1;
        param_general.flag_cirrus
    case "old_local_test"
        param_general.image_name = 'cubeW28';
        spectral_downsampling = 20;
        spatial_downsampling = 4;
        coverage_path = "data/vla_7.95h_dt10s.uvw256.mat";
        param_model.ncores_data = 2;
        param_model.Qx = 1;
        param_model.Qy = 1;
        param_general.flag_cirrus
    otherwise
        error("Unknown experiment type");
end

reference_cube_path = fullfile(data_path, strcat(param_general.image_name, '.fits'));
info        = fitsinfo(reference_cube_path);
rowend      = info.PrimaryData.Size(1);
colend      = info.PrimaryData.Size(2);
sliceend    = info.PrimaryData.Size(3);
if strcmp(param_model.algo_version, 'sara')
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
n_channels = floor(sliceend / spectral_downsampling);

[Ny, Nx, nchans] = size(x0);
N = Nx * Ny;
X0 = reshape(x0, [N, nchans]);

% frequency used to generate the reference cubes
nu0 = 2.052e9;  % starting freq
dnu = 16e6;  % freq step
L = 100;  % number of channels
nu_vect = [nu0 (dnu * (1:L - 1) + nu0)];
frequencies = nu_vect(1:floor(L / n_channels):end);

clear reference_cube_path info rowend colend sliceend;
clear spatial_downsampling spectral_downsampling;

disp('MNRAS configuration');
disp(['Algorithm version: ', param_model.algo_version]);
disp(['Reference image: ', param_general.image_name]);
% disp(['nchannels: ', num2str(n_channels)]);
disp(['Number of facets param_model.Qy x param_model.Qx : ', num2str(param_model.Qy), ' x ', num2str(param_model.Qx)]);
disp(['Number of spectral facets param_model.Qc : ', num2str(param_model.Qc)]);
disp(['Overlap fraction: ', strjoin(strsplit(num2str(param_model.overlap_fraction)), ', ')]);

%% Auxiliary function needed to select the appropriate workers
% (only needed for 'hs' and 'fhs' algorithms)
switch param_model.algo_version
    case 'sara'
        data_worker_id = @(k) k;
    case 'hs'
        data_worker_id = @(k) k;
    case 'fhs'
        data_worker_id = @(k) k + param_model.Qx * param_model.Qy;
    otherwise
        error('Undefined param_model.algo_version');
end

%% Get faceting parameter (spectral + spatial)
% fix faceting parameters in case these are not consistent with the
% selected algorithm
if strcmp(param_model.algo_version, 'sara')
    param_model.window_type = 'none';
    param_model.Qc = n_channels;
    param_model.Qx = 1;
    param_model.Qy = 1;
elseif strcmp(param_model.algo_version, 'hs')
    param_model.window_type = 'none';
    param_model.Qx = 1;
    param_model.Qy = 1;
end
Q = param_model.Qx * param_model.Qy;

% convert fraction of overlap between consecutive facets into a number of pixels
overlap_size = get_overlap_size([Ny, Nx], [param_model.Qy, param_model.Qx], param_model.overlap_fraction);
disp(['Number of pixels in overlap: ', strjoin(strsplit(num2str(overlap_size)), ' x ')]);

% index of the spectral channels involved in the subcube
interleaved_channels = split_range_interleaved(param_model.Qc, n_channels);
subcube_channels = interleaved_channels{ind};

% index of channels from the subcube to be handled on each data worker
rg_c = split_range(param_model.ncores_data, nchans);

% frequencies associated with the current subcube
fc = frequencies(subcube_channels);
fmax = frequencies(end);

% extract ground truth subcube
if param_model.Qc > 1 && ind > 0 && ~strcmp(param_model.algo_version, 'sara')
    x0 = x0(:, :, subcube_channels);
    nchans = size(x0, 3);
    X0 = reshape(x0, Nx * Ny, nchans);
end

%% Setup name of results file
data_name_function = @(nchannels) strcat('y_', ...
    param_synth.exp_type, '_', param_general.image_name, '_srf=', num2str(param_synth.superresolution_factor), ...
    '_Ny=', num2str(Ny), '_Nx=', num2str(Nx), '_L=', num2str(nchannels), '_snr=', num2str(param_synth.isnr), ...
    '.mat');

temp_results_name = @(nchannels) strcat(param_synth.exp_type, '_', param_general.image_name, '_', ...
    param_model.algo_version, '_', param_model.window_type, '_srf=', num2str(param_synth.superresolution_factor), ...
    '_Ny=', num2str(Ny), '_Nx=', num2str(Nx), '_L=', num2str(nchannels), ...
    '_Qy=', num2str(param_model.Qy), '_Qx=', num2str(param_model.Qx), '_Qc=', num2str(param_model.Qc), ...
    '_ind=', num2str(ind), '_g=', num2str(param_model.gamma), '_gb=', num2str(param_model.gamma_bar), ...
    '_overlap=', strjoin(strsplit(num2str(param_model.overlap_fraction)), '_'));

warm_start = @(nchannels) strcat(temp_results_name(nchannels), '_rw=', num2str(param_model.warmstart_iteration), '.mat');

data_name = data_name_function(n_channels);

%% Define problem configuration (rng, nufft, preconditioning, blocking,
% NNLS (epsilon estimation), SARA dictionary)
% parameters_problem;

%% Generate/load uv-coverage
% generating u-v coverage
% ! reminder uv-coverage and weighting
% https://casa.nrao.edu/Release4.1.0/doc/UserMan/UserMansu259.html
if param_synth.flag_generateCoverage
    cov_type = 'vlaa';
    p = 0.5;
    dl = 1.1;
    hrs = 5;
    na = 27; % for vlaa
    M = na * (na - 1) / 2;
    % Fixing Mt = 0.5 N, take T = 0.5 N / M : na = 27 for vla
    T = floor(p * (Nx * Ny) / M); % should be > 1
    [u, v, ~] = generate_uv_coverage(T, hrs, dl, cov_type);
    u = u(:) * fc(1) / fmax;
    v = v(:) * fc(1) / fmax;
    fitswrite([u, v, ones(numel(u), 1)], coverage_path);
    disp(coverage_path);
else
    disp(strcat("Loading coverage: ", coverage_path));

    % VLA configuration
    % A. 762775 -> 3
    % B. 268448 -> 2
    % C. 202957 -> 1
    % D. 47750 -> 0

    if strcmp(param_synth.exp_type, "spectral")
        load(coverage_path, 'uvw', 'obsId');
        size(uvw);
        u1 = uvw(obsId == 2, 1) * fmax / speed_of_light;
        v1 = uvw(obsId == 2, 2) * fmax / speed_of_light;
        clear obsId;
    else
        % ! normalize u,v coverage w.r.t. the highest frequency (i.e., uv expressed in
        % units of the smallest wavelength, associated with the highest frequency)
        load(coverage_path, 'uvw');
        size(uvw);
        u1 = uvw(:, 1) * fmax / speed_of_light;
        v1 = uvw(:, 2) * fmax / speed_of_light;
%         load(coverage_path, 'uvw', 'obsId');
%         size(uvw)
%         u1 = uvw(obsId==3, 1)*fmax/speed_of_light;
%         v1 = uvw(obsId==3, 2)*fmax/speed_of_light;
%         clear obsId
    end
    bmax = max(sqrt(u1.^2 + v1.^2));

    % cellsize = 3600*180/(param_synth.superresolution_factor*2*pi*bmax); % in arcsec
    u = u1 * pi / (param_synth.superresolution_factor * bmax);
    v = v1 * pi / (param_synth.superresolution_factor * bmax);
    size(u);
    disp('Coverage loaded successfully');
    clear uvw u1 v1;
end

%% Setup parpool
cirrus_cluster = util_set_parpool(param_model.algo_version, param_model.ncores_data, param_model.Qx * param_model.Qy, param_general.flag_cirrus);

%% Setup measurement operator
switch param_model.algo_version
    case 'sara'
        if param_model.flagDR
            % ! define Sigma (weight matrix involved in DR)
            % ! define G as the holographic matrix
        else
            [A, At, G, W, aW] = util_gen_measurement_operator(u, v, ...
                param_precond, param_blocking, fc, fmax, Nx, Ny, ...
                param_nufft.Kx, param_nufft.Ky, param_nufft.ox, param_nufft.oy);
            Sigma = [];
        end

        % if ~param_model.flagDR
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
        if strcmp(param_model.algo_version, 'hs')
            spmd
                local_fc = fc(rg_c(labindex, 1):rg_c(labindex, 2));
                if param_model.flagDR
                    % ! define Sigma (weight matrix involved in DR)
                    % ! define G as the holographic matrix
                else
                    % ! ideally, simplify irt nufft interface to do so
                    [A, At, G, W, aW] = util_gen_measurement_operator(u, v, ...
                    param_precond, param_blocking, local_fc, fmax, Nx, Ny, param_nufft.Kx, param_nufft.Ky, param_nufft.ox, param_nufft.oy, param_nufft.kernel);
                    Sigma = [];
                end
            end
        else
            spmd
                % define operator on data workers only
                if labindex > Q
                    local_fc = fc(rg_c(labindex - Q, 1):rg_c(labindex - Q, 2));
                    if param_model.flagDR
                        % ! define Sigma (weight matrix involved in DR)
                        % ! define G as the holographic matrix
                    else
                        % ! ideally, simplify irt nufft interface to do so
                        [A, At, G, W, aW] = util_gen_measurement_operator(u, v, ...
                        param_precond, param_blocking, local_fc, fmax, Nx, Ny, param_nufft.Kx, param_nufft.Ky, param_nufft.ox, param_nufft.oy, param_nufft.kernel);
                        Sigma = [];
                    end
                end
            end
        end
        clear local_fc;
end

%% Free memory
clear param_blocking param_precond;

%% Generate/load visibilities (generate only full spectral dataset)
% only generatr data in 'hs' or 'fhs' configuration (otherwise, load the data)
datafile = matfile(fullfile(results_path, data_name));

switch param_model.algo_version
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

        for k = 1:param_model.ncores_data
            y{data_worker_id(k)} = datafile.y(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1);
            epsilons{data_worker_id(k)} = datafile.epsilons(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1);
            sigma_noise{data_worker_id(k)} = datafile.sigma_noise(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1);
        end
        global_sigma_noise = datafile.sigma_noise;
end
disp('Data loaded successfully');

%% Compute operator norm
if strcmp(param_model.algo_version, 'sara')
    if param_general.flag_compute_operator_norm
        [Anorm, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-8, 200);

        save(fullfile(results_path, ...
            strcat('Anorm_', param_model.algo_version, ...
            '_Ny=', num2str(Ny), '_Nx=', num2str(Nx), ...
            '_L=', num2str(n_channels), ...
            '_Qc=', num2str(param_model.Qc), '_ind=', num2str(ind), ...
            '_ch=', num2str(ind), '.mat')), ...
            '-v7.3', 'Anorm', 'squared_operator_norm', 'rel_var', ...
            'squared_operator_norm_precond', 'rel_var_precond');
        clear rel_var;
    else
        load(fullfile(results_path, ...
            strcat('Anorm_', param_model.algo_version, ...
            '_Ny=', num2str(Ny), '_Nx=', num2str(Nx), ...
            '_L=', num2str(n_channels), ...
            '_Qc=', num2str(param_model.Qc), '_ind=', num2str(ind), ...
            '_ch=', num2str(ind), '.mat')), ...
            'Anorm', 'squared_operator_norm_precond', 'squared_operator_norm');
    end
else
    if param_general.flag_compute_operator_norm
        spmd
            if labindex > param_model.Qx * param_model.Qy * strcmp(param_model.algo_version, 'fhs')
                [An, squared_operator_norm, rel_var, squared_operator_norm_precond, rel_var_precond] = util_operator_norm(G, W, A, At, aW, Ny, Nx, 1e-8, 200);
            end
        end

        % save operator norm from the different subcubes into a single .mat
        % file
        opnormfile = matfile(fullfile(results_path, strcat('Anorm', ...
            '_Ny=', num2str(Ny), '_Nx=', num2str(Nx), ...
            '_L=', num2str(n_channels), '.mat')), 'Writable', true);

        opnormfile.squared_operator_norm = zeros(nchans, 1);
        opnormfile.rel_var = zeros(nchans, 1);
        opnormfile.squared_operator_norm_precond = zeros(nchans, 1);
        opnormfile.rel_var_precond = zeros(nchans, 1);

        Anorm = 0;
        for k = 1:param_model.ncores_data
            opnormfile.squared_operator_norm(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1) = squared_operator_norm{data_worker_id(k)};
            opnormfile.rel_var(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1) = rel_var{data_worker_id(k)};

            opnormfile.squared_operator_norm_precond(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1) = squared_operator_norm_precond{data_worker_id(k)};
            opnormfile.rel_var_precond(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1) = rel_var_precond{data_worker_id(k)};

            Anorm = max(Anorm, An{data_worker_id(k)});
        end
        squared_operator_norm = opnormfile.squared_operator_norm;
        squared_operator_norm_precond = opnormfile.squared_operator_norm_precond;
        clear An rel_var rel_var_precond;

    else
        opnormfile = matfile(fullfile(results_path, strcat('Anorm', ...
            '_Ny=', num2str(Ny), '_Nx=', num2str(Nx), ...
            '_L=', num2str(n_channels), '.mat')));

        squared_operator_norm_precond = opnormfile.squared_operator_norm_precond(subcube_channels, 1);
        rel_var_precond = opnormfile.rel_var_precond(subcube_channels, 1);
        Anorm = max(squared_operator_norm_precond .* (1 + rel_var_precond));
        squared_operator_norm = opnormfile.squared_operator_norm(subcube_channels, 1);

        % squared_operator_norm = Composite();
        % for k = 1:param_model.ncores_data
        %     squared_operator_norm{data_worker_id(k)} = opnormfile.squared_operator_norm(subcube_channels(rg_c(k, 1)):subcube_channels(rg_c(k, 2)), 1);
        % end
    end
end

fprintf('Convergence parameter (measurement operator): %e \n', Anorm);

%% Regularization parameters and solver

% estimate noise level (set regularization parameters to the same value)
% compute sig and sig_bar (estimate of the "noise level" in "SVD" and
% SARA space) involved in the reweighting scheme

if strcmp(param_model.algo_version, 'sara')

    % SARA dicionary (created out of the solver for SARA)
    dwtmode('zpd');
    [Psi1, Psit1] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
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
    mu = param_model.gamma * sig;
    fprintf('Noise level: sig = %e\n', sig);
    fprintf('Additional multiplicative regularisation factor param_model.gamma = %e\n', param_model.gamma);
    fprintf('Regularization parameter mu = %e\n', mu);
    fprintf('Algo: %s, alpha = %.4e, mu = %.4e, sig = %.4e\n', param_model.algo_version, param_model.gamma, mu, sig);
end

if strcmp(param_model.algo_version, 'hs') || strcmp(param_model.algo_version, 'fhs')

    % noise level / regularization parameter
    [sig, sig_bar, mu_chi, sig_chi, sig_sara] = ...
        compute_noise_level(Ny, Nx, nchans, global_sigma_noise, ...
        param_model.algo_version, param_model.Qx, param_model.Qy, overlap_size, squared_operator_norm);

    % apply multiplicative factor for the regularization parameters (if needed)
    mu_bar = param_model.gamma_bar * sig_bar;
    mu = param_model.gamma * sig;
    fprintf('mu_chi = %.4e, sig_chi = %.4e, sig_sara = %.4e\n', mu_chi, sig_chi, sig_sara);
    fprintf('Noise levels: sig = %.4e, sig_bar = [%.4e, %.4e]\n', sig, min(sig_bar), max(sig_bar));
    fprintf('Additional multiplicative actors param_model.gamma = %.4e, param_model.gamma_bar = %.4e\n', param_model.gamma, param_model.gamma_bar);
    fprintf('Regularization parameters: mu = %.4e, mu_bar = %.4e\n', mu, mu_bar);
    fprintf('Algo: %s, param_model.gamma = %.4e, param_model.gamma_bar = %.4e, mu = %.4e, mu_bar = [%.4e, %.4e]\n', param_model.algo_version, param_model.gamma, param_model.gamma_bar, mu, min(mu_bar), max(mu_bar));
end

% * general
% parameters_solver;
% estimate of the noise level in SARA space
param_solver.reweighting_sig = sig;
if ~strcmp(param_model.algo_version, 'sara')
    % estimate of the noise level in "SVD" spaces
    param_solver.reweighting_sig_bar = sig_bar;
end
% bound on the norm of the Identity operator
param_solver.nu0 = 1;
% bound on the norm of the operator Psi
param_solver.nu1 = 1;
% upper bound on the norm of the measurement operator
param_solver.nu2 = squared_operator_norm_precond;
% regularization parameter nuclear norm
if ~strcmp(param_model.algo_version, 'sara')
    param_solver.gamma0 = mu_bar;
end
% regularization parameter l21-norm (soft th parameter)
% ! for SARA, take the value given as an input to the solver
param_solver.gamma = mu;
% id of the cube to be reconstructed (if spectral faceting is active)
param_solver.cube_id = ind;

%%
% TODO: update solver interface
name_checkpoint = fullfile(auxiliary_path, temp_results_name(n_channels));
name_warmstart = fullfile(auxiliary_path, warm_start(n_channels));

if param_general.flag_solve_minimization
    %%
    if strcmp(param_model.algo_version, 'sara')
        disp('SARA');
        disp('-----------------------------------------');

        % ! in this case, param_model.ncores_data corresponds
        % ! to the number of workers for the wavelet transform (9 maximum)
        xsol = sara(y, epsilons, A, At, aW, G, W, Psi, Psit, ...
                    param_solver, name_warmstart, name_checkpoint, param_model.gamma, ...
                    param_model.flagDR, Sigma, [], x0);

        mkdir('results/');

        fitswrite(xsol, fullfile(auxiliary_path, strcat('x_', param_general.image_name, '_', param_model.algo_version, ...
        '_srf=', num2str(param_synth.superresolution_factor), ...
        '_', param_model.window_type, ...
        '_Qy=', num2str(param_model.Qy), '_Qx=', num2str(param_model.Qx), '_Qc=', num2str(param_model.Qc), ...
        '_ind=', num2str(ind), ...
        '_gam=', num2str(param_model.gamma), ...
        '.fits')));
    else
        %%

        % spectral tesselation (non-overlapping)
        % ! to be updated tonight (need to be careful about the different variables needed + implicit parallelization conventions)
        cell_c_chunks = cell(param_model.ncores_data, 1); % ! to check
        for k = 1:param_model.ncores_data
            cell_c_chunks{k} = rg_c(k, 1):rg_c(k, 2);
        end
        
        % ! SARA dictionary prior (hard-coded)
        % depth of the wavelet decompositions
        nlevel = 4;
        % wavelet dictionaries
        % ! always specify Dirac basis ('self') in last position if ever used
        wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};
        % length of the filters (by convention, 0 corresponds to the Dirac basis)
        filter_length = [2 * (1:8)'; 0];

        %%
        switch param_model.algo_version
            case 'hs'
                disp('HyperSARA');
                disp('-----------------------------------------');
                xsol = hyperSARA(y, epsilons, ...
                    A, At, aW, G, W, param_solver, ...
                    param_model.ncores_data, wlt_basis, nlevel, cell_c_chunks, ...
                    nchans, Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, param_model.flagDR, Sigma, ...
                    [], X0);
            case 'fhs'
                disp('Faceted HyperSARA');
                disp('-----------------------------------------');
                xsol = facetHyperSARA(y, epsilons, ...
                    A, At, aW, G, W, param_solver, param_model.Qx, param_model.Qy, param_model.ncores_data, ...
                    wlt_basis, filter_length, nlevel, param_model.window_type, ...
                    cell_c_chunks, nchans, overlap_size, param_model.gamma, param_model.gamma_bar, ...
                    Ny, Nx, param_nufft.oy, param_nufft.ox, ...
                    name_warmstart, name_checkpoint, param_model.flagDR, Sigma, ...
                    [], X0);
            otherwise
                error('Unknown solver version.');
        end

        mkdir('results/');
        fitswrite(xsol, fullfile(auxiliary_path, strcat('x_', param_general.image_name, '_', param_model.algo_version, ...
            '_', param_model.window_type, ...
            '_srf=', num2str(param_synth.superresolution_factor), ...
            '_Qy=', num2str(param_model.Qy), '_Qx=', num2str(param_model.Qx), '_Qc=', num2str(param_model.Qc), ...
            '_ind=', num2str(ind), ...
            '_gam=', num2str(param_model.gamma), '_gambar=', num2str(param_model.gamma_bar), ...
            '_overlap=', strjoin(strsplit(num2str(param_model.overlap_fraction)), '_'), ...
            '.fits')));
    end
end
