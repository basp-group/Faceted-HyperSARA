function nnls_dr_reduced_data(chInd, reduction_version, realdatablocks, enable_klargestpercent, fouRed_gamma)
% chInd = 1;
% reduction_version = 2;
fprintf('Channel number %d\n', chInd);
fprintf('Reduction version %d\n', reduction_version);
fprintf('Data blocks: %d\n', realdatablocks);

addpath ../lib/utils/
addpath ../fouRed/
addpath ../lib/operators
addpath ../lib/nufft
% addpath ../data_mnras_dr

param_real_data.image_size_Nx = 2560;
param_real_data.image_size_Ny = 1536;
% nSpw = 16;          % number of spectral channels per MS file
% nChannels = 1*nSpw; % total number of "virtual" channels (i.e., after
% % concatenation) for the real dataset considered
nBlocks = realdatablocks;        % number of data blocks (needs to be known beforehand,
% quite restrictive here), change l.70 accordingly
klargestpercent = 20;
FT2 = @(x) fftshift(fft2(ifftshift(x))) / sqrt(numel(x));

%% Config parameters
Nx = param_real_data.image_size_Nx;
Ny = param_real_data.image_size_Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

% %% Preconditioning parameters
% param_precond.N = N;       % number of pixels in the image
% param_precond.Nox = ox*Nx; % number of pixels in the image
% param_precond.Noy = oy*Ny; % number of pixels in the image
% param_precond.gen_uniform_weight_matrix = 1; % weighting type
% param_precond.uniform_weight_sub_pixels = 1;

%% Fourier reduction parameters
param_fouRed.enable_klargestpercent = enable_klargestpercent;
param_fouRed.klargestpercent = klargestpercent;
param_fouRed.enable_estimatethreshold = ~enable_klargestpercent;
param_fouRed.gamma = fouRed_gamma;             % By using threshold estimation, the optimal theshold reads as gamma * sigma / ||x||_2
param_fouRed.diagthresholdepsilon = 0;
param_fouRed.fastCov = 1;

fprintf('Enable k-largest percentage: %d\n', param_fouRed.enable_klargestpercent);
fprintf('Enable automatic threshold: %d\n', param_fouRed.enable_estimatethreshold);
%% parameter NNLS
param_nnls.verbose = 2;       % print log or not
param_nnls.rel_obj = 1e-4;    % stopping criterion
param_nnls.max_iter = 200;     % max number of iterations 1000
param_nnls.sol_steps = [inf]; % saves images at the given iterations
param_nnls.beta = 1;

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

param_block_structure.use_manual_partitioning = 0;

if param_block_structure.use_manual_partitioning == 1
    param_block.size = 90000; % to be changed (make sure nBlocks in the end...)
    param_block.snapshot = 0;
end

% data files on cirrus
% new_file_y = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_y.mat');
% new_file_u = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_u.mat');
% new_file_v = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_v.mat');
% new_file_nW = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_nW.mat');
% data files on workstation
% new_file_y = matfile('/home/basphw/mjiang/Data/mjiang/extract_real_data/CYG_2b_y.mat');
% new_file_u = matfile('/home/basphw/mjiang/Data/mjiang/extract_real_data/CYG_2b_u.mat');
% new_file_v = matfile('/home/basphw/mjiang/Data/mjiang/extract_real_data/CYG_2b_v.mat');
% new_file_nW = matfile('/home/basphw/mjiang/Data/mjiang/extract_real_data/CYG_2b_nW.mat');
% data files in EPFL
% new_file_y = matfile('/Users/ming/workspace/Git/extract_real_data/CYG_y.mat');
% new_file_u = matfile('/Users/ming/workspace/Git/extract_real_data/CYG_u.mat');
% new_file_v = matfile('/Users/ming/workspace/Git/extract_real_data/CYG_v.mat');
% new_file_nW = matfile('/Users/ming/workspace/Git/extract_real_data/CYG_nW.mat');
% load('/Users/ming/workspace/Git/extract_real_data/CYG_data_raw_ind=6.mat', 'y', 'uw', 'vw', 'nWw')
load(['/lustre/home/shared/sc004/cyg_data_sub/CYG_data_cal_', num2str(realdatablocks), 'b_ch32_ind=6.mat'], 'yb', 'G')
load(['/lustre/home/shared/sc004/cyg_data_sub/res_', num2str(realdatablocks), 'b_ch32_ind=6.mat'], 'res')
%% Reduction 
epsilon = cell(1, 1);

% define operators
% parpool(6)
[A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);

% instantiate variables for DR
H = cell(1, 1);  % holographic matrices
yT = cell(1, 1); % reduced data
T = cell(1, 1);  % preconditioning matrix (inverse of singular values)
Wm = cell(1, 1); % mask DR
aW = cell(1, 1); % preconditioner

% solve NNLS per block / estimate epsilon / reduce data
% y_tmp = new_file_y.y(chInd,1);
% u_tmp = new_file_u.u(chInd,1);
% v_tmp = new_file_v.v(chInd,1);
% nW_tmp = new_file_nW.nW(chInd,1);
% y_tmp = y_tmp{1};
% u_tmp = u_tmp{1};
% v_tmp = v_tmp{1};
% nW_tmp = nW_tmp{1};

y_tmp = yb{1};
% u_tmp = uw;
% v_tmp = vw;
% nW_tmp = nWw;
G_tmp = G{1};
res_tmp = res{1};

H{1} = cell(numel(y_tmp), 1);
T{1} = cell(numel(y_tmp), 1);
yTl = cell(numel(y_tmp), 1);
Wm{1} = cell(numel(y_tmp), 1);
aWl = cell(numel(y_tmp), 1);
norm_res = cell(numel(y_tmp), 1);

Hl = H{1};
Wml = Wm{1};
Tl = T{1};

if param_fouRed.enable_estimatethreshold
    if ~isfield(param_fouRed, 'gamma') 
        param_fouRed.gamma = 3; 
    end
    fprintf('Threshold level: %d sigma\n', param_fouRed.gamma);
    p = normcdf([-param_fouRed.gamma param_fouRed.gamma]);
    prob = p(2) - p(1);
end

diagthresholdepsilon = param_fouRed.diagthresholdepsilon;

% parpool(nBlocks)
for j = 1:nBlocks
    
%     load G_H_2560_1536.mat
    
%     [A, At, G, ~] = op_p_nufft([v_tmp(j) u_tmp(j)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW_tmp(j));

    Hl{j} = (G_tmp{j}')*G_tmp{j}; % progressively write to disk? (possibly huge...)
    
    if reduction_version == 1    
    %     fast matrix probing (using psf)
        dirac2D = zeros(Ny, Nx);
        dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;
        PSF = operatorIpsf(dirac2D, A, At, Hl{j}, [oy*Ny, ox*Nx]);
        covariancemat = FT2(PSF);
        d_mat = abs((covariancemat(:)));
    elseif reduction_version == 2
        d_mat = full(abs(diag(Hl{j})));
    end
    
    if param_fouRed.enable_klargestpercent
        th = 0;
        Mask = (d_mat > th);
%         Mask = (d_mat >= prctile(d_mat,100-param_fouRed.klargestpercent));
    elseif param_fouRed.enable_estimatethreshold
    %     estimate threshold
        d_mat_sort = sort(d_mat);
        d_mat_sort_cumsum = cumsum(d_mat_sort);
        d_mat_sum = d_mat_sort_cumsum(end); % whole energy of d_mat
        th_energy = d_mat_sum * (1 - prob); % threshold according to the k-sigma rule
        th = d_mat_sort(find(d_mat_sort_cumsum >= th_energy, 1, 'first'));
        Mask = (d_mat >= th);
    end
    
    percentage = sum(Mask(:)) / numel(Mask) * 100;
    fprintf('Block %d of channel %d, %f%% of data are kept, threshold=%f\n', j, chInd, percentage, th);
    
    d_mat = d_mat(Mask);

    Tl{j} = max(diagthresholdepsilon, d_mat);  % ensures that inverting the values will not explode in computation
    Tl{j} = 1./sqrt(Tl{j});
    Wml{j} = Mask;
    
    if reduction_version == 1
%         aWl{j} = 1./Tl{j};
        aWl{j} = 1;
        yTl{j} = dataReduce(y_tmp{j}, G_tmp{j}', At, Tl{j}, Wml{j});
    elseif reduction_version == 2
        aWl{j} = 1;
%         yTl{j} = dataReduce_degrid(y_tmp{j}, G{1}', Tl{j}, Wml{j});
%         yTl{j} = dataReduce_degrid(y_tmp{j}, G_tmp{j}', Tl{j}, Wml{j});
        resRed = dataReduce_degrid(res_tmp{j}, G_tmp{j}', Tl{j}, Wml{j});       % !!! for calibrated data
        norm_res{j} = norm(resRed);
    end
    
%     [~, norm_res{j}] = fb_dr_nnls(yTl{j}, A, At, Hl{j}, Tl{j}, Wml{j}, param_nnls, reduction_version);
    fprintf('Block %d of channel %d, estimated epsilon: %f\n',j, chInd, norm_res{j})

end
% T{1} = Tl;
% Wm{1} = Wml;
% aW{1} = aWl;
% yT{1} = yTl;
% H{1} = Hl;
epsilon{1} = norm_res;

% % save on workstation
% save(['/home/basphw/mjiang/Data/mjiang/real_data_dr/CYG_old_epsilon=', num2str(chInd), '.mat'],'-v7.3', 'epsilon');
% save(['/home/basphw/mjiang/Data/mjiang/real_data_dr/CYG_old_yT=', num2str(chInd), '.mat'],'-v7.3', 'yT');
% save(['/home/basphw/mjiang/Data/mjiang/real_data_dr/CYG_old_DR=', num2str(chInd), '.mat'],'-v7.3', 'H', 'T', 'aW', 'Wm');

% % save in EPFL
% save(['./CYG_epsilon=', num2str(chInd), '.mat'],'-v7.3', 'epsilon');
% save(['./CYG_H_9b=', num2str(chInd), '.mat'],'-v7.3', 'H');
% save(['./CYG_DR=', num2str(chInd), '.mat'],'-v7.3', 'H', 'T', 'aW', 'Wm');
% save(['./CYG_DR_9b=', num2str(chInd), '.mat'],'-v7.3', 'yT', 'T', 'aW', 'Wm', 'epsilon');

% save on cirrus
% Hfilename = ['/lustre/home/shared/sc004/dr_', num2str(realdatablocks), 'b_result_real_data/CYG_H_cal_', num2str(realdatablocks), 'b_ind6=', num2str(chInd), '.mat'];
% if ~isfile(Hfilename)
%     save(Hfilename, '-v7.3', 'H')
% end
% DRfilename = ['/lustre/home/shared/sc004/dr_', num2str(realdatablocks), 'b_result_real_data/CYG_DR_cal_', num2str(realdatablocks), 'b_ind6_fouRed',...
%     num2str(reduction_version), '_th', num2str(fouRed_gamma),'=', num2str(chInd), '.mat'];
% if ~isfile(DRfilename)
%     save(DRfilename, '-v7.3', 'yT', 'T', 'aW', 'Wm', 'epsilon');
% end
Epsfilename = ['/lustre/home/shared/sc004/dr_', num2str(realdatablocks), 'b_result_real_data/CYG_eps_cal_', num2str(realdatablocks), 'b_ind6_fouRed',...
    num2str(reduction_version), '_th', num2str(fouRed_gamma),'=', num2str(chInd), '.mat'];
if ~isfile(Epsfilename)
    save(Epsfilename, '-v7.3', 'epsilon');
end
fprintf('Dimensionality reduction and epsilon estimation are finished\n')
% save(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_G_epsilon=', num2str(chInd), '.mat'],'-v7.3', 'epsilon_new');
% save(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_G_yT=', num2str(chInd), '.mat'],'-v7.3', 'yT');
% save(['/lustre/home/shared/sc004/dr_2b_result_real_data/CYG_G_DR=', num2str(chInd), '.mat'],'-v7.3', 'T', 'aW', 'Wm');
