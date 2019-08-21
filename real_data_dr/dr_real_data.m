function dr_real_data(chInd)
fprintf('Channel number %d', chInd);

addpath ../lib/utils/
addpath ../fouRed/
addpath ../lib/operators
addpath ../lib/nufft
% addpath ../data_mnras_dr

param_real_data.image_size_Nx = 2560;
param_real_data.image_size_Ny = 1536;
nSpw = 16;          % number of spectral channels per MS file
nChannels = 2*nSpw; % total number of "virtual" channels (i.e., after
% concatenation) for the real dataset considered
nBlocks = 2;        % number of data blocks (needs to be known beforehand,
% quite restrictive here), change l.70 accordingly
klargestpercent = 20;
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
param_precond.Nox = ox*Nx; % number of pixels in the image
param_precond.Noy = oy*Ny; % number of pixels in the image
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
param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi/4), 0, icdf('norm', 0.75, 0, pi/4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)

param_block_structure.use_manual_partitioning = 1;

if param_block_structure.use_manual_partitioning == 1
    param_block.size = 90000; % to be changed (make sure nBlocks in the end...)
    param_block.snapshot = 0;
end

% data files on cirrus
new_file_y = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_y.mat');
new_file_u = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_u.mat');
new_file_v = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_v.mat');
new_file_nW = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_nW.mat');
new_file_res = matfile(['/lustre/home/shared/sc004/nnls_real_data/CYG_2b_res_nnls=', num2str(chInd),'.mat']);
% data files on workstation
% new_file_y = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_y.mat');
% new_file_u = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_u.mat');
% new_file_v = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_v.mat');
% new_file_nW = matfile('/home/basphw/mjiang/Data/pthouvenin/basp_sharing/Ming/data_mnras_dr/CYG_nW.mat');
% new_file_res = matfile(['/home/basphw/mjiang/Data/mjiang/real_data_dr/CYG_res_nnls=', num2str(chInd),'.mat'])

%% Estimate epsilon with NNLS on each data block

%     if generate_eps_nnls
%         % set up res_nnls to be saved on disk
%         res_nnls = cell(nChannels, 1);
%         for l = 1:numel(res_nnls)
%             res_nnls{l} = cell(nBlocks, 1);
%         end
%         save('CYG_res_nnls.mat', 'res_nnls', '-v7.3');
%         clear res_nnls
%         
%         % parameter NNLS
%         param_nnls.verbose = 2;       % print log or not
%         param_nnls.rel_obj = 1e-5;    % stopping criterion
%         param_nnls.max_iter = 1000;     % max number of iterations 1000
%         param_nnls.sol_steps = [inf]; % saves images at the given iterations
%         param_nnls.beta = 1;
%     end
epsilon = cell(1, 1);

% define operators
% parpool(6)
[A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);

% instantiate variables for DR
H = cell(1, 1);  % holographic matrices
yT = cell(1, 1); % reduced data
T = cell(1, 1);  % preconditioning matrix (inverse of singular values)
Wm = cell(1, 1); % mask DR
W = cell(1, 1);  % masks for the Fourier plane corresponding to
% data blocks
aW = cell(1, 1); % preconditioner

% solve NNLS per block / estimate epsilon / reduce data
y_tmp = new_file_y.y(chInd,1);
u_tmp = new_file_u.u(chInd,1);
v_tmp = new_file_v.v(chInd,1);
nW_tmp = new_file_nW.nW(chInd,1);
res_tmp = new_file_res.res_nnls(1,1);
y_tmp = y_tmp{1};
u_tmp = u_tmp{1};
v_tmp = v_tmp{1};
nW_tmp = nW_tmp{1};
res_tmp = res_tmp{1};

H{1} = cell(numel(u_tmp), 1);
T{1} = cell(numel(u_tmp), 1);
yTl = cell(numel(u_tmp), 1);
Wm{1} = cell(numel(u_tmp), 1);
aWl = cell(numel(u_tmp), 1);
eps_ = cell(numel(u_tmp), 1);

Hl = H{1};
Wml = Wm{1};
Tl = T{1};

gamma = param_fouRed.gamma;
diagthresholdepsilon = param_fouRed.diagthresholdepsilon;

parpool(nBlocks)
parfor j = 1:nBlocks
    [~, ~, G, ~] = op_p_nufft([v_tmp(j) u_tmp(j)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW_tmp(j));

%     eps(j, l) = norm(res_tmp{j});
    Hl{j} = (G{1}')*G{1}; % progressively write to disk? (possibly huge...)

%     estimate threshold
    dirty2 = norm(operatorPhit(y_tmp{j}, G{1}', At), 'fro');

%     fast matrix probing (using psf)
    dirac2D = zeros(Ny, Nx);
    dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;
    PSF = operatorIpsf(dirac2D, A, At, Hl{j}, [oy*Ny, ox*Nx]);
    covariancemat = FT2(PSF);
    d_mat = abs(real(covariancemat(:)));

    rn = FT2(At(G{1}'*res_tmp{j}));
    th_dirty = gamma * std(rn(:)) / (dirty2 / max(d_mat(:)));
    Mask = (d_mat >= th_dirty);
    d_mat = d_mat(Mask);

    Tl{j} = max(diagthresholdepsilon, d_mat);  % ensures that inverting the values will not explode in computation
    Tl{j} = 1./sqrt(Tl{j});
    Wml{j} = Mask;
    aWl{j} = 1./Tl{j};

%     reduce the data block and residual
    yTl{j} = dataReduce(y_tmp{j}, G{1}', At, Tl{j}, Wml{j});
    reduced_res = dataReduce(res_tmp{j}, G{1}', At, Tl{j}, Wml{j});
    eps_{j} = norm(reduced_res(:), 2);
end
T{1} = Tl;
Wm{1} = Wml;
aW{1} = aWl;
yT{1} = yTl;
H{1} = Hl;
epsilon{1} = eps_;

% % save on workstation
% save(['/home/basphw/mjiang/Data/mjiang/real_data_dr/CYG_epsilon=', num2str(chInd), '.mat'],'-v7.3', 'epsilon');
% save(['/home/basphw/mjiang/Data/mjiang/real_data_dr/CYG_yT=', num2str(chInd), '.mat'],'-v7.3', 'yT');
% save(['/home/basphw/mjiang/Data/mjiang/real_data_dr/CYG_DR=', num2str(chInd), '.mat'],'-v7.3', 'H', 'T', 'aW', 'Wm');

% save on cirrus
save(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_epsilon=', num2str(chInd), '.mat'],'-v7.3', 'epsilon');
save(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_yT=', num2str(chInd), '.mat'],'-v7.3', 'yT');
save(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_DR=', num2str(chInd), '.mat'],'-v7.3', 'H', 'T', 'aW', 'Wm');
