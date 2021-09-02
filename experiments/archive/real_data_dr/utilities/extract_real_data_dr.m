% script to extract necessary data from raw data
%% Real data extraction
addpath ../lib/utils/;
addpath ../fouRed/;
addpath ../lib/operators;
addpath ../lib/measurement-operator/nufft;
% addpath ../data_mnras_dr

visibility_file_name = '/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_data_raw_ind=';
param_real_data.image_size_Nx = 2560;
param_real_data.image_size_Ny = 1536;
nSpw = 16;          % number of spectral channels per MS file
nChannels = 2 * nSpw; % total number of "virtual" channels (i.e., after
% concatenation) for the real dataset considered
nBlocks = 9;        % number of data blocks (needs to be known beforehand,
% quite restrictive here), change l.70 accordingly
klargestpercent = 20;
extract_real_data = false;
generate_eps_nnls = true;
reduce_data = true;
FT2 = @(x) fftshift(fft2(ifftshift(x)));

%% Config parameters
Nx = param_real_data.image_size_Nx;
Ny = param_real_data.image_size_Ny;
N = Nx * Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

%% Preconditioning parameters
param_precond.N = N;       % number of pixels in the image
param_precond.Nox = ox * Nx; % number of pixels in the image
param_precond.Noy = oy * Ny; % number of pixels in the image
param_precond.gen_uniform_weight_matrix = 1; % weighting type
param_precond.uniform_weight_sub_pixels = 1;

%% Fourier reduction parameters
param_fouRed.enable_klargestpercent = 1;
param_fouRed.klargestpercent = klargestpercent;
param_fouRed.enable_estimatethreshold = 0;
param_fouRed.gamma = 3;             % By using threshold estimation, the optimal theshold reads as gamma * sigma / ||x||_2
param_fouRed.diagthresholdepsilon = 1e-10;
param_fouRed.covmatfileexists = 0;
param_fouRed.covmatfile = 'covariancemat.mat';
param_fouRed.fastCov = 1;

%% Block structure
regenerate_block_structure = 1;

param_block_structure.use_density_partitioning = 0;
param_block_structure.density_partitioning_no = 1;

param_block_structure.use_uniform_partitioning = 0;
param_block_structure.uniform_partitioning_no = 4;

param_block_structure.use_equal_partitioning = 0;
param_block_structure.equal_partitioning_no = 1;

param_block_structure.use_manual_frequency_partitioning = 0;
% sparam.fpartition = [pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [0, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.25*pi, 0, 0.25*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-64/256*pi, 0, 64/256*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi / 4), 0, icdf('norm', 0.75, 0, pi / 4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)

param_block_structure.use_manual_partitioning = 1;

if param_block_structure.use_manual_partitioning == 1
    param_block.size = 90000; % to be changed (make sure nBlocks in the end...)
    param_block.snapshot = 0;
end

%% Load/extract real data from band-interleaved .mat files
if extract_real_data
    % visibilities
    y = cell(nChannels, 1);
    for l = 1:numel(y)
        y{l} = cell(nBlocks, 1);
    end
    save('CYG_y.mat', 'y', '-v7.3');
    clear y;

    % u
    u = cell(nChannels, 1);
    for l = 1:numel(u)
        u{l} = cell(nBlocks, 1);
    end
    save('CYG_u.mat', 'u', '-v7.3');
    clear u;

    % v
    v = cell(nChannels, 1);
    for l = 1:numel(v)
        v{l} = cell(nBlocks, 1);
    end
    save('CYG_v.mat', 'v', '-v7.3');
    clear v;

    % nWw: scaling NUFFT
    nW = cell(nChannels, 1);
    for l = 1:numel(nW)
        nW{l} = cell(nBlocks, 1);
    end
    save('CYG_nW.mat', 'nW', '-v7.3');
    clear nW;

    %     % position of data for different configs inside each visibility vector
    %     % (e.g., A, B, C ...)
    %     pos = cell(nChannels, 1);
    %     for l = 1:numel(pos)
    %         pos{l} = cell(nBlocks, 1);
    %     end
    %     save('CYG_pos.mat', 'pos', '-v7.3');
    %     clear pos;

    %     % acquisition time
    %     time = cell(nChannels, 1);
    %     for l = 1:numel(time)
    %         time{l} = cell(nBlocks, 1);
    %     end
    %     save('CYG_time.mat', 'time', '-v7.3');
    %     clear time;

    new_file_y = matfile('CYG_y.mat', 'Writable', true);
    new_file_u = matfile('CYG_u.mat', 'Writable', true);
    new_file_v = matfile('CYG_v.mat', 'Writable', true);
    new_file_nW = matfile('CYG_nW.mat', 'Writable', true);

    for l = 1:nChannels
        % concatenate the visibility frequencies and the associated u/v
        % points (make sure same number of blocks)
        u_tmp = new_file_u.u(l, 1);
        v_tmp = new_file_v.v(l, 1);
        y_tmp = new_file_y.y(l, 1);
        nW_tmp = new_file_nW.nW(l, 1);
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);
            y = cell2mat(file.y(1, l));
            pos = cell2mat(file.pos(1, l));
            time = cell2mat(file.time_(1, l));
            uw = cell2mat(file.uw(1, l));
            vw = cell2mat(file.vw(1, l));
            nWw = cell2mat(file.nWw(1, l));

            % blocking: set the blocks structure (per frequency)
            param_block.pos = pos;
            out_block = util_time_based_block_sp_ar(uw, time, param_block);
            param_block_structure.partition = out_block.partition;
            aWw = util_gen_preconditioning_matrix(uw, vw, param_precond);
            [u1, v1, ~, uvidx1, aW, nW1] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);

            for m = 1:numel(u1)
                u_tmp{1}{m} = [u_tmp{1}{m}; u1{m}];
                v_tmp{1}{m} = [v_tmp{1}{m}; v1{m}];
                y_tmp{1}{m} = [y_tmp{1}{m}; y(uvidx1{m})];
                nW_tmp{1}{m} = [nW_tmp{1}{m}; nW1{m}];
            end
        end
        new_file_u.u(l, 1) = u_tmp;
        new_file_v.v(l, 1) = v_tmp;
        new_file_y.y(l, 1) = y_tmp;
        new_file_nW.nW(l, 1) = nW_tmp;
    end
else
    new_file_y = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_y.mat');
    new_file_u = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_u.mat');
    new_file_v = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_v.mat');
    new_file_nW = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_nW.mat');
end

%% Estimate epsilon with NNLS on each data block
if reduce_data

    if generate_eps_nnls
        % set up res_nnls to be saved on disk
        res_nnls = cell(nChannels, 1);
        for l = 1:numel(res_nnls)
            res_nnls{l} = cell(nBlocks, 1);
        end
        save('CYG_res_nnls.mat', 'res_nnls', '-v7.3');
        clear res_nnls;

        % parameter NNLS
        param_nnls.verbose = 2;       % print log or not
        param_nnls.rel_obj = 1e-5;    % stopping criterion
        param_nnls.max_iter = 1000;     % max number of iterations 1000
        param_nnls.sol_steps = [inf]; % saves images at the given iterations
        param_nnls.beta = 1;
    end
    new_file_res = matfile('CYG_res_nnls.mat', 'Writable', true);
    epsilon = cell(nChannels, 1);
    eps = zeros(nBlocks, nChannels);

    % define operators
    parpool(6);
    [A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy * Ny ox * Nx], [Ny / 2 Nx / 2]);

    % instantiate variables for DR
    H = cell(nChannels, 1);  % holographic matrices
    yT = cell(nChannels, 1); % reduced data
    T = cell(nChannels, 1);  % preconditioning matrix (inverse of singular values)
    Wm = cell(nChannels, 1); % mask DR
    W = cell(nChannels, 1);  % masks for the Fourier plane corresponding to
    % data blocks
    aW = cell(nChannels, 1); % preconditioner

    % solve NNLS per block / estimate epsilon / reduce data
    for l = 1:nChannels
        y_tmp = new_file_y.y(l, 1);
        u_tmp = new_file_u.u(l, 1);
        v_tmp = new_file_v.v(l, 1);
        nW_tmp = new_file_nW.nW(l, 1);
        res_tmp = new_file_res.res_nnls(l, 1);
        y_tmp = y_tmp{1};
        u_tmp = u_tmp{1};
        v_tmp = v_tmp{1};
        nW_tmp = nW_tmp{1};
        res_tmp = res_tmp{1};

        H{l} = cell(numel(u_tmp), 1);
        T{l} = cell(numel(u_tmp), 1);
        yTl = cell(numel(u_tmp), 1);
        Wm{l} = cell(numel(u_tmp), 1);
        aWl = cell(numel(u_tmp), 1);
        eps_ = cell(numel(u_tmp), 1);

        Hl = H{l};
        Wml = Wm{l};
        Tl = T{l};

        for j = 1:nBlocks
            [~, ~, G, W1] = op_p_nufft([v_tmp(j) u_tmp(j)], [Ny Nx], [Ky Kx], [oy * Ny ox * Nx], [Ny / 2 Nx / 2], nW_tmp(j));
            if generate_eps_nnls
                [sol, ~] = fb_nnls_blocks(y_tmp{j}, A, At, G{1}(:, W1{1}), W1{1}, param_nnls); % check format of "sol" (matrix)
                res_tmp{j} = y_tmp{j} - G{1} * A(sol); % check computation of the residual
            end
%             eps(j, l) = norm(res_tmp{j});
%             Hl{j} = (G{1}')*G{1}; % progressively write to disk? (possibly huge...)
%
%             estimate threshold
%             dirty2 = norm(operatorPhit(y_tmp{j}, G{1}', At) / sqrt(N));
%
%             fast matrix probing (using psf)
%             dirac2D = zeros(Ny, Nx);
%             dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;
%             PSF = operatorIpsf(dirac2D, A, At, Hl{j}, [oy*Ny, ox*Nx]);
%             covariancemat = FT2(PSF);
%             d_mat = abs(real(covariancemat(:)));
%
%             rn = FT2(At(G{1}'*res_tmp{j}));
%             th_dirty = param_fouRed.gamma * std(rn(:)) / (dirty2 / max(d_mat(:)));
%             Mask = (d_mat >= th_dirty);
%             d_mat = d_mat(Mask);
%
%             Tl{j} = max(param_fouRed.diagthresholdepsilon, d_mat);  % ensures that inverting the values will not explode in computation
%             Tl{j} = 1./sqrt(Tl{j});
%             Wml{j} = Mask;
%             aWl{j} = 1./Tl{j};
%
%             reduce the data block and residual
%             yTl{j} = dataReduce(y_tmp{j}, G{1}', At, Tl{j}, Wml{j});
%             reduced_res = dataReduce(res_tmp{j}, G{1}', At, Tl{j}, Wml{j});
%             eps_{j} = norm(reduced_res(:), 2);
        end
%         T{l} = Tl;
%         Wm{l} = Wml;
%         aW{l} = aWl;
%         yT{l} = yTl;
%         H{l} = Hl;
%         epsilon{l} = eps_;
        new_file_res.res_nnls(l, 1) = {res_tmp};
    end
%     save('CYG_epsilon.mat','-v7.3', 'epsilon', 'eps');
%     save('CYG_yT.mat','-v7.3', 'yT');
%     save('CYG_DR.mat','-v7.3', 'H', 'T', 'aW', 'Wm');
else
    load('CYG_DR.mat');
    load('CYG_yT.mat', '-v7.3');
    load('CYG_epsilon.mat');
end
