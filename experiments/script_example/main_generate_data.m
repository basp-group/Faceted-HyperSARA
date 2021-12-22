function main_generate_data(json_filename, ind)

%%
% TODO: update util_gen_measurement_operator to enable kaiser kernels
% json_filename = "setup_matlab.json";
% ind = 1;
format compact;

% seed for the random number generation (used only when generating data)
seed = 0;

% Read .json configuration file
[speed_of_light, param_general, param_model, ~, param_nufft, ...
 param_blocking, param_precond, param_synth, ~] = ...
    read_json_configuration(json_filename, ind);

% %% Debug parameters
% param_general.image_name = 'W28_512'; %'cygASband_Cube_H'; %'W28_512';
% param_synth.exp_type = 'local_test'; % 'spectral', 'spatial', 'test'
% param_synth.flag_generateCoverage = 0;
% param_model.flagDR = 0;
% 
% param_model.ncores_data = 2; % number of cores assigned to the data fidelity terms (groups of channels)
% param_general.path_uv_coverage = "data/vla_7.95h_dt10s.uvw256.mat" ;%"data/msSpecs.mat"; % "data/vla_7.95h_dt10s.uvw256.mat";
% param_synth.isnr = 50;
% 
% param_synth.superresolution_factor = 2;
% param_general.flag_cirrus = false;
% param_nufft.kernel = 'minmax:tuned'; % 'kaiser' (for real data), 'minmax:tuned'

addpath ../../lib/operators/;
addpath ../../lib/measurement-operator/nufft/;
addpath ../../lib/measurement-operator/lib/operators/;
addpath ../../lib/generate_data;
addpath ../../lib/measurement-operator/lib/utils/;
addpath ../../lib/faceted-wavelet-transform/src;
addpath ../../lib/utils/;
addpath ../../data/;
addpath ../../src/;

% setting paths to results and reference image cube
data_path = '../../data/';
results_path = fullfile('results', strcat(param_general.image_name, '_', param_synth.exp_type));
mkdir(data_path);
mkdir(results_path);

%%
% ! load the appropriate portion of the reference image cube
% ! if spectral faceting is active, just load the interesting portion of 
% ! the full image cube
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
        param_general.path_uv_coverage = "data/vla_7.95h_dt10s.uvw256.mat";
        param_model.ncores_data = 2;
        param_general.flag_cirrus
    case "old_local_test"
        param_general.image_name = 'cubeW28';
        spectral_downsampling = 20;
        spatial_downsampling = 4;
        param_general.path_uv_coverage = "data/vla_7.95h_dt10s.uvw256.mat";
        param_model.ncores_data = 2;
        param_general.flag_cirrus
    otherwise
        error("Unknown experiment type");
end

reference_cube_path = fullfile(data_path, strcat(param_general.image_name, '.fits'));
info        = fitsinfo(reference_cube_path);
rowend      = info.PrimaryData.Size(1);
colend      = info.PrimaryData.Size(2);
sliceend    = info.PrimaryData.Size(3);

x0 = fitsread(reference_cube_path, 'primary', ...
            'Info', info, ...
            'PixelRegion', {[1 spatial_downsampling rowend], ...
            [1 spatial_downsampling colend], ...
            [1 spectral_downsampling sliceend]});

n_channels = floor(sliceend / spectral_downsampling);

[Ny, Nx, nchans] = size(x0);
input_snr = param_synth.isnr * ones(nchans, 1);  % input SNR (in dB)

disp(['Number of channels considered: ', num2str(nchans)]);

% frequency used to generate the reference cubes
nu0 = 2.052e9;  % starting freq
dnu = 16e6;  % freq step
L = 100;  % number of channels
nu_vect = [nu0 (dnu * (1:L - 1) + nu0)];
frequencies = nu_vect(1:floor(L / n_channels):end);

clear reference_cube_path info rowend colend sliceend;
clear spatial_downsampling spectral_downsampling;

disp('Synthetic data generation');
disp(['Reference image: ', param_general.image_name]);
disp(['Input SNR: ', num2str(param_synth.isnr)]);

%% Auxiliary function needed to select the appropriate workers
% index of channels from the subcube to be handled on each data worker
rg_c = split_range(param_model.ncores_data, nchans);

% frequencies associated with the current subcube
fmax = frequencies(end);

%% Setup name of data file
data_name_function = @(nchannels) strcat('y_', ...
    param_synth.exp_type, '_', param_general.image_name, '_srf=', num2str(param_synth.superresolution_factor), ...
    '_Ny=', num2str(Ny), '_Nx=', num2str(Nx), '_L=', num2str(nchannels), ...
    '_snr=', num2str(param_synth.isnr), '.mat');

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
    u = u(:) * frequencies(1) / fmax;
    v = v(:) * frequencies(1) / fmax;
    fitswrite([u, v, ones(numel(u), 1)], param_general.path_uv_coverage);
    disp(param_general.path_uv_coverage);
else
    disp(strcat("Loading coverage: ", param_general.path_uv_coverage));
    
    % VLA configuration
    % A. 762775 -> 3
    % B. 268448 -> 2
    % C. 202957 -> 1
    % D. 47750 -> 0

    if strcmp(param_synth.exp_type, "spectral")
        load(param_general.path_uv_coverage, 'uvw', 'obsId');
        size(uvw);
        u1 = uvw(obsId == 2, 1) * fmax / speed_of_light;
        v1 = uvw(obsId == 2, 2) * fmax / speed_of_light;
        clear obsId;
    else
        % ! normalize u,v coverage w.r.t. the highest frequency (i.e., uv expressed in
        % units of the smallest wavelenght, associated with the highest frequency)
        load(param_general.path_uv_coverage, 'uvw');
        size(uvw);
        u1 = uvw(:, 1) * fmax / speed_of_light;
        v1 = uvw(:, 2) * fmax / speed_of_light;
%         load(param_general.path_uv_coverage, 'uvw', 'obsId');
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
cirrus_cluster = util_set_parpool('hs', param_model.ncores_data, 1, param_general.flag_cirrus);

%% Setup measurement operator
% create the measurement operator operator in parallel (depending on
% the algorithm used)
spmd
    local_fc = frequencies(rg_c(labindex, 1):rg_c(labindex, 2));
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
clear local_fc;

%% Free memory
clear param_blocking param_precond;

%% Generate/load visibilities (generate only full spectral dataset)
param_l2_ball.type = 'sigma';
param_l2_ball.sigma_ball = 2;

% https://coderedirect.com/questions/287901/random-number-generator-matlab-with-multiple-cpus
% https://fr.mathworks.com/help/matlab/math/creating-and-controlling-a-random-number-stream.html?searchHighlight=random%20number%20streams&s_tid=srchtitle#brvku_2
% may not be strictly equivalent to the data initially obtained for the
% mnras paper
[rng_stream{1:param_model.ncores_data}] = ...
    RandStream.create('threefry4x64_20', 'Seed', seed, 'NumStreams', ...
    param_model.ncores_data);

spmd
    [y0, y, Ml, ~, sigma_noise, ~] = util_gen_measurements_snr( ...
        x0(:, :, rg_c(labindex, 1):rg_c(labindex, 2)), G, W, A, ...
        input_snr(rg_c(labindex, 1):rg_c(labindex, 2)), rng_stream{labindex});
    [~, epsilons] = util_gen_data_fidelity_bounds2(y, Ml, .../
        param_l2_ball, sigma_noise);
end

% save parameters (matfile solution)
datafile = matfile(fullfile(results_path, data_name), ...
    'Writable', true);
datafile.y0 = cell(nchans, 1);
datafile.y = cell(nchans, 1);
datafile.epsilons = cell(nchans, 1);
datafile.sigma_noise = zeros(nchans, 1);

% ! need to convert from local to global channel index (i.e., index 
% ! into the full spectral dimension by using 
% ! subcube_channels(rg_c(k, 1))
for k = 1:param_model.ncores_data
    datafile.y0(rg_c(k, 1):rg_c(k, 2), 1) = y0{k};
    datafile.y(rg_c(k, 1):rg_c(k, 2), 1) = y{k};
    datafile.epsilons(rg_c(k, 1):rg_c(k, 2), 1) = epsilons{k};
    datafile.sigma_noise(rg_c(k, 1):rg_c(k, 2), 1) = sigma_noise{k};
end