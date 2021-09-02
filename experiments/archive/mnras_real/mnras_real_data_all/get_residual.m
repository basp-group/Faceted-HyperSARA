% function get_residual(mat_filenames, fits_filenames, extract_raw_data, ...
% normalizeResiduals)
% Compute residual image (if missing from a given run)
%    mat_filenames: cell contaning the name of the files to be addressed
%    fits_filenames: name of the corresponding fits files
%
%%

% ! to be adapted with the full path to the files, with extension
% ! to be commented of script used as a function
% extract_raw_data = false;
% normalizeResiduals = false; % option to normalize the residual images or not
% results_path = fullfile('results/');
% base_mat_filenames = ...
% base_mat_filenames = ...
% mat_filenames = cell(nfiles, 1);
% fits_filenames = cell(nfiles, 1);
% for k = 1:nfiles
%     mat_filenames{k} = ...
%     fits_filenames{k} = ...
% end
% !----

format compact;
% raise an error if the two cells do not have the same number of elements
if abs(numel(mat_filenames) - numel(fits_filenames)) > 0
    error("get_residual:incorrectSize", "Error: the input cells `mat_filenames` and `fits_filenames` should have the same size.");
end
addpath ../lib/operators/;
addpath ../lib/measurement-operator/nufft/;
addpath ../lib/utils/;
addpath ../lib/faceted-wavelet-transform/src;
addpath ../data/;
addpath ../src_mnras/;

%% Default config parameters (can be kept as is)
% gridding parameters
Nx = 2560;
Ny = 1536;
N = Nx * Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 7; % number of neighbours for nufft
Ky = 7; % number of neighbours for nufft

% preconditioning parameters
param_precond.N = N;       % number of pixels in the image
param_precond.Nox = ox * Nx; % number of Fourier points (oversampled plane)
param_precond.Noy = oy * Ny;
param_precond.gen_uniform_weight_matrix = 1; % set weighting type
param_precond.uniform_weight_sub_pixels = 1;

% block structure
param_block_structure.use_density_partitioning = 0;
param_block_structure.density_partitioning_no = 1;
param_block_structure.use_uniform_partitioning = 0;
param_block_structure.uniform_partitioning_no = 4;
param_block_structure.use_equal_partitioning = 0;
param_block_structure.equal_partitioning_no = 1;
param_block_structure.use_manual_frequency_partitioning = 0;
param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi / 4), 0, icdf('norm', 0.75, 0, pi / 4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_block_structure.use_manual_partitioning = 1;

%% Extracting data
% ! @Abdullah: is it still needed, or did you saved the data in a .mat file?
% ! by default the following should be fine
visibility_file_name = {'5', '7'};
configuration_file_name = {'C', 'A'};

% load real_data;
if extract_raw_data
    [y, uw1, vw1, nWw, f, time_, pos] = util_load_real_vis_data(visibility_file_name, configuration_file_name, ind);
    [uw, vw, pixel_size] = scale_uv(visibility_file_name, configuration_file_name, uw1, vw1);
    save(['CYG_data_raw_ind=' num2str(ind) '.mat'], '-v7.3', 'y', 'uw', 'vw', 'nWw', 'f', 'time_', 'pos', 'pixel_size');
else
    load(['CYG_data_raw_ind=' num2str(ind) '.mat']);
end

ch = [1:17, 19, 21:32];
ii = 1;

for i = ch
    i;

    %% compute weights
    [aWw] = util_gen_preconditioning_matrix(uw{i}, vw{i}, param_precond);

    % set the blocks structure
    param_block_structure.partition = pos{i};
    [u, v, ~, uvidx, aW{ii}, nW] = util_gen_block_structure(uw{i}, vw{i}, aWw, nWw{i}, param_block_structure);
    u;
    v;
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    tstart = tic;

    % compute A & At with Kx = Ky = 7
    [A, At, ~, ~] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy * Ny ox * Nx], [Ny / 2 Nx / 2], nW); % Kx = Ky = 7

    % load the G matrices and the estimated vis.
    if i < 17
        load(['/lustre/home/shared/sc004/AD_Update_G_Codes/results/SubCube' num2str(ind) '/WB5-' num2str(i) '/PreProcStruct.mat']);
    else
        load(['/lustre/home/shared/sc004/AD_Update_G_Codes/results/SubCube' num2str(ind) '/WB7-' num2str(i) '/PreProcStruct.mat']);
    end
    PreProcStruct.Gw;
    Gw = [PreProcStruct.Gw{2}; PreProcStruct.Gw{1}]; % concatenate the configurations C first and then A
    % param_HSI.xsol_Arwa(:,:,ii) = PreProcStruct.SolInit;
    clear PreProcStruct;
    disp('Gw');
    size(Gw);

    G{ii} = cell(length(u), 1);
    yb{ii} = cell(length(u), 1);
    W{ii}  = cell(length(u), 1); %%%%%%%% ARWA
    for j = 1:length(u)
        G{ii}{j} = Gw(uvidx{j}, :);
        W{ii}{j} = G{ii}{j}' * ones(size(G{ii}{j}, 1), 1) ~= 0;
        G{ii}{j} = G{ii}{j}(:, W{ii}{j});
        yb{ii}{j} = y{i}(uvidx{j});
    end
    clear Gw res_est_full;

    ii = ii + 1;
end

%% Free memory
clear y u v uv_mat uv_mat1 uvidx uw vw ant1 ant2 aWw nW nWw out_block;
clear param_block_structure param_precond;
clear u v u1 v1 uw vw aWw nW nWw param_block_structure param_precond;
nChannels = size(yb, 2);

for k = 1:numel(mat_filenames)
    disp(mat_filenames{k});
    % load image from backup file
    mfile = matfile(mat_filenames{k}, 'Writable', true);
    xsol = mfile.xsol;
    % compute residual image
    res = zeros(Ny, Nx, nChannels);
    dirac_2d = zeros(Ny, Nx);
    dirac_2d(floor(Ny / 2) + 1, floor(Nx / 2) + 1) = 1;
    F_dirac = A(dirac_2d);
    psf_peak = zeros(nChannels);

    for l = 1:nChannels
        z = zeros(oy * Ny * ox * Nx, 1);
        Fx = A(xsol(:, :, l));
        z_psf = zeros(oy * Ny * ox * Nx, 1);

        for b = 1:numel(yb{l})
            z(W{l}{b}) = z(W{l}{b}) + G{l}{b}' * (yb{l}{b} - (G{l}{b} * Fx(W{l}{b})));
            z_psf(W{l}{b}) = z_psf(W{l}{b}) + G{l}{b}' * (G{l}{b} * F_dirac(W{l}{b}));
        end

        res(:, :, l) = real(At(z));
        psf = real(At(z_psf));
        psf_peak(l) = max(psf(:));
        if normalizeResiduals
            res(:, :, l) = res(:, :, l) / psf_peak(l);
        end

    end
    % update residual image in .mat backup file
    mfile.res = res;
    mfile.psf_peak = psf_peak;
    if normalizeResiduals
        m.normalizedResiduals = false;
    else
        m.normalizedResiduals = true;
    end
    % save residual image as .fits file
    fitswrite(res, fits_filenames{k});
end
