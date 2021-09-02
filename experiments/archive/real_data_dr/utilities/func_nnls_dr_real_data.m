function func_nnls_dr_real_data(chInd, reduction_version, realdatablocks, enable_klargestpercent, fouRed_gamma)
% Generate and save non-calibrated DR data from raw data
% -------------------------------------------------------------------------%
% Input:
% > chInd: vector of channel indices [L]
% > reduction_version: reduction version,
%       1: old one (not used any more)
%       2: current one
% > realdatablocks: number of data blocks
% > enable_klargestpercent: activate strategy of reduction of removing
% smallest k percent singular values, otherwise strategy of reduction of
% removing (statitical notion) k-sigma smallest singular values
% > fouRed_gamma: level of reduction

fprintf('Channel number %d\n', chInd);
fprintf('Reduction version %d\n', reduction_version);
fprintf('Data blocks: %d\n', realdatablocks);

addpath ../lib/utils/;
addpath ../fouRed/;
addpath ../lib/operators;
addpath ../lib/measurement-operator/nufft/;
% addpath ../data_mnras_dr

Nx = 2560;
Ny = 1536;
nBlocks = realdatablocks;        % number of data blocks for real data
% klargestpercent = 20;
FT2 = @(x) fftshift(fft2(ifftshift(x))) / sqrt(numel(x));

%% Config parameters
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 7; % number of neighbours for nufft
Ky = 7; % number of neighbours for nufft

% %% Preconditioning parameters
% param_precond.N = N;       % number of pixels in the image
% param_precond.Nox = ox*Nx; % number of pixels in the image
% param_precond.Noy = oy*Ny; % number of pixels in the image
% param_precond.gen_uniform_weight_matrix = 1; % weighting type
% param_precond.uniform_weight_sub_pixels = 1;

%% Fourier reduction parameters
param_fouRed.enable_klargestpercent = enable_klargestpercent;
param_fouRed.enable_estimatethreshold = ~enable_klargestpercent;
param_fouRed.gamma = fouRed_gamma;  % reduction parameter
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

%% data files on cirrus
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
% load(['/Users/ming/workspace/Git/extract_real_data/CYG_data_cal_', num2str(realdatablocks), 'b_ch', num2str(chInd),'_ind=6.mat'], 'yb', 'G', 'W')
% load(['/Users/ming/workspace/Git/extract_real_data/res_', num2str(realdatablocks), 'b_ch', num2str(chInd),'_ind=6.mat'], 'res')
load(['/lustre/home/shared/sc004/cyg_data_sub/CYG_data_cal_', num2str(realdatablocks), 'b_ch', num2str(chInd), '_ind=6.mat'], 'yb', 'G', 'W');
load(['/lustre/home/shared/sc004/cyg_data_sub/res_', num2str(realdatablocks), 'b_ch', num2str(chInd), '_ind=6.mat'], 'res');
%% Reduction
epsilon = cell(1, 1);

% define operators
% parpool(6)
[A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy * Ny ox * Nx], [Ny / 2 Nx / 2]);

% instantiate variables for DR
H = cell(1, 1);  % holographic matrices
yT = cell(1, 1); % reduced data
T = cell(1, 1);  % normalization (inverse of singular values)
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
W_tmp = W{1};
res_tmp = res{1};       % precomputed residuals for calibrated data

H{1} = cell(numel(y_tmp), 1);
yT{1} = cell(numel(y_tmp), 1);
T{1} = cell(numel(y_tmp), 1);
Wm{1} = cell(numel(y_tmp), 1);
aW{1} = cell(numel(y_tmp), 1);
norm_res = cell(numel(y_tmp), 1);
% sing{1} = cell(numel(y_tmp), 1);

Hl = H{1};
yTl = yT{1};
Tl = T{1};
Wml = Wm{1};
aWl = aW{1};

if param_fouRed.enable_estimatethreshold
    if ~isfield(param_fouRed, 'gamma')
        param_fouRed.gamma = 30;
    end
    fprintf('Threshold level: remove %d percentile\n', param_fouRed.gamma);
%     p = normcdf([-param_fouRed.gamma param_fouRed.gamma]);
%     prob = p(2) - p(1);
    prob = 1 - param_fouRed.gamma / 100;
end

precision = 1e-16;

nb_raw = 0;
nb_red = 0;
for j = 1:nBlocks

%     [A, At, G, ~] = op_p_nufft([v_tmp(j) u_tmp(j)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW_tmp(j));

    Hl{j} = (G_tmp{j}') * G_tmp{j}; % progressively write to disk? (possibly huge...)
    peak = max(max(abs(Hl{j})));
    Hl{j} = Hl{j} .* (abs(Hl{j}) > peak * precision);

    if reduction_version == 1
    %     fast matrix probing (using psf)
        dirac2D = zeros(Ny, Nx);
        dirac2D(ceil((Ny + 1) / 2), ceil((Nx + 1) / 2)) = 1;
        PSF = operatorIpsf(dirac2D, A, At, Hl{j}, [oy * Ny, ox * Nx]);
        covariancemat = FT2(PSF);
        d_mat = abs(covariancemat(:));
    elseif reduction_version == 2
        d_mat = full(abs(diag(Hl{j})));
    end

    if param_fouRed.enable_klargestpercent
        Mask = (d_mat > param_fouRed.gamma);
%         Mask = (d_mat >= prctile(d_mat,100-param_fouRed.klargestpercent));
    elseif param_fouRed.enable_estimatethreshold
    %     estimate threshold
        d_mat_sort = sort(d_mat);
        d_mat_sort_cumsum = cumsum(d_mat_sort);
        d_mat_sum = d_mat_sort_cumsum(end); % whole energy of d_mat
        th_energy = d_mat_sum * (1 - prob); % threshold according to the k-sigma rule
        th = d_mat_sort(find(d_mat_sort_cumsum >= th_energy, 1, 'first'));
        Mask = (d_mat >= th);
        ind_nz = d_mat > 0;       % non-zero singular values
    end

    kept = sum(Mask(:));
    total = numel(Mask);
    percentage = kept / total * 100;
    nb_red = nb_red + kept;
    nb_raw = nb_raw + numel(y_tmp{j});
    fprintf('Block %d of channel %d: %d non-zero singular values, %d over %d, or %f%% of data are kept, threshold=%f\n', j, chInd, sum(ind_nz), kept, total, percentage, th);

    d_mat = d_mat(Mask);
    Hl{j} = Hl{j}(Mask, :);

    Tl{j} = d_mat;
%     Tl{j} = max(diagthresholdepsilon, d_mat);  % ensures that inverting the values will not explode in computation
    Tl{j} = 1 ./ sqrt(Tl{j});
    Wml{j} = Mask;

    if reduction_version == 1
        aWl{j} = 1;
        yTl{j} = dataReduce(y_tmp{j}, G_tmp{j}', W_tmp{j}, At, Tl{j}, Wml{j});
        resRed = dataReduce(res_tmp{j}, G_tmp{j}', W_tmp{j}, At, Tl{j}, Wml{j});       % !!! for calibrated data, residuals are known
        norm_res{j} = norm(resRed);
    elseif reduction_version == 2
        aWl{j} = 1;
        Gt = G_tmp{j}';
        yTl{j} = Tl{j} .* (Gt(Wml{j}, :) * y_tmp{j});
        resRed = Tl{j} .* (Gt(Wml{j}, :) * res_tmp{j});        % !!! for calibrated data, residuals are known
        norm_res{j} = norm(resRed);
    end

%     [~, norm_res{j}] = fb_dr_nnls(yTl{j}, A, At, Hl{j}, W_tmp{j}, Tl{j}, Wml{j}, param_nnls, reduction_version);
    fprintf('Block %d of channel %d, estimated epsilon: %f\n', j, chInd, norm_res{j});

end
fprintf('Raw data size: %d, reduced data size: %d\n', nb_raw, nb_red);
fprintf('H memory: \n');
whos Hl;

H{1} = Hl;
yT{1} = yTl;
T{1} = Tl;
aW{1} = aWl;
Wm{1} = Wml;
epsilon{1} = norm_res;

DRfilename = ['/lustre/home/shared/sc004/dr_', num2str(realdatablocks), 'b_result_real_data/CYG_DR_cal_', num2str(realdatablocks), 'b_ind6_fouRed', ...
    num2str(reduction_version), '_perc', num2str(fouRed_gamma), '=', num2str(chInd), '.mat'];
if ~isfile(DRfilename)
    save(DRfilename, '-v7.3', 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'epsilon');
end
fprintf('Dimensionality reduction and epsilon estimation are finished\n');
