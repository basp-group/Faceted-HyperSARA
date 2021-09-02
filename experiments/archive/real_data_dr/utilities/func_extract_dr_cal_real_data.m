function func_extract_dr_cal_real_data(chInd, subInd, reduction_version, realdatablocks, enable_klargestpercent, fouRed_gamma, fouRed_type)
% Generate and save calibrated DR data from raw data and calibrated G matrices
% -------------------------------------------------------------------------%
% Input:
% > chInd: vector of channel indices [L]
% > subInd: vector of interlaced subproblem indices [S]
% > reduction_version: reduction version,
%       1: old one (not used any more)
%       2: current one
% > realdatablocks: number of data blocks
% > enable_klargestpercent: activate strategy of reduction of removing
% smallest k percent singular values, otherwise strategy of reduction of
% removing (statitical notion) k-sigma smallest singular values
% > fouRed_gamma: level of reduction
% > fouRed_type: type of reduction, supplementary information for
% fouRed_gamma
%       1: remove smallest singular values based on "fouRed_gamma" percentage
%       2: remove smallest singular values based on "fouRed_gamma"-sigma

addpath ../../lib/utils/;
addpath ../../fouRed/;
addpath ../../lib/operators;
addpath ../../lib/measurement-operator/nufft;

fprintf('Channel number: %d\n', chInd);
fprintf('Index number: %d\n', subInd);
fprintf('Reduction version %d\n', reduction_version);
if fouRed_type == 1
    typeStr = 'perc';
elseif fouRed_type == 2
    typeStr = 'th';
end
fprintf('Data blocks: %d\n', realdatablocks);

%% Data directory
if chInd == 31          % ch31 is reflagged
    datadir = '/lustre/home/shared/sc004/SubCubes_ch31/';
else
    datadir = '/lustre/home/shared/sc004/AD_Data_Subcubes/SubCubes/';
end
gmatdir = '/lustre/home/shared/sc004/AD_Update_G_Codes/results/';

%% real data aggregation
Nx = 2560;
Ny = 1536;

xsol = zeros(Ny, Nx);

Gw = cell(realdatablocks, 1);
Wl = cell(realdatablocks, 1);
res = cell(realdatablocks, 1);
yb = cell(realdatablocks, 1);

%%%%%% Trick to use more RAM on cirrus %%%%%%
% delete(gcp('nocreate'));
% numworkers = 2;
% cirrus_cluster = parcluster('cirrus R2019a');
% cirrus_cluster.NumWorkers = numworkers;
% cirrus_cluster.NumThreads = 1;
% ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
% if cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads > ncores
%     exit(1);
% end
% parpool(cirrus_cluster, numworkers);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:length(subInd)
    fprintf('\nIndex number: %d\n', subInd(j));
    if chInd < 17
        kerl = 5;
    else
        kerl = 7;
    end
    % Calibrated G matrices
    gmatfile = [gmatdir, 'SubCube', num2str(subInd(j)), '/WB', num2str(kerl), '-', num2str(chInd), '/PreProcStruct.mat'];
    fprintf('Read G matrix file: %s\n', gmatfile);
    tmp = load(gmatfile);
    xsol = xsol + tmp.PreProcStruct.SolInit;
    for k = 1:realdatablocks
        Gw{k} = [Gw{k}; tmp.PreProcStruct.Gw{k}];
        res{k} = [res{k}; tmp.PreProcStruct.ResiudalModelDataInit{k}];
    end

    % Raw data
    if chInd == 31
        datafile = [datadir, 'CYG', num2str(j), '_ch31.mat'];
    else
        datafile = [datadir, 'CYG', num2str(j), '.mat'];
    end
    fprintf('Read data file: %s\n', datafile);
    tmp = load(datafile);
    for k = 1:realdatablocks
        if chInd == 31
            yb{k} = [yb{k}; tmp.y_I{1}{k}(tmp.y_I{1}{k} ~= 0)];
        else
            yb{k} = [yb{k}; tmp.y_I{chInd}{k}(tmp.y_I{chInd}{k} ~= 0)];
        end
    end
end
xsol = xsol / length(subInd);   % average of initial xsol

if realdatablocks == 2
    nb_a = numel(yb{1});
    nb_c = numel(yb{2});
end
nb_t = nb_a + nb_c;
fprintf('\nNumber of A-config data points: %d, number of C-config data points: %d, number of total data points: %d\n', nb_a, nb_c, nb_t);

for k = 1:realdatablocks
    Wl{k} = Gw{k}' * ones(size(Gw{k}, 1), 1) ~= 0;       % remove zero columns of G matrices
    Gw{k} = Gw{k}(:, Wl{k});
end
fprintf('G matrix memory:\n');
whos Gw;

clear tmp;

%% Reduction
FT2 = @(x) fftshift(fft2(ifftshift(x))) / sqrt(numel(x));
% Config parameters
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 7; % number of neighbours for nufft
Ky = 7; % number of neighbours for nufft

% Fourier reduction parameters
param_fouRed.enable_klargestpercent = enable_klargestpercent;
param_fouRed.enable_estimatethreshold = ~enable_klargestpercent;
param_fouRed.gamma = fouRed_gamma;  % reduction parameter
param_fouRed.diagthresholdepsilon = 0;
param_fouRed.fastCov = 1;

fprintf('Dimensionality reduction begins ...\n');
fprintf('Enable k-largest percentage: %d\n', param_fouRed.enable_klargestpercent);
fprintf('Enable automatic threshold: %d\n', param_fouRed.enable_estimatethreshold);

% define operators
% parpool(6)
[A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy * Ny ox * Nx], [Ny / 2 Nx / 2]);

% instantiate variables for DR
H = cell(1, 1);  % holographic matrices
yT = cell(1, 1); % reduced data
T = cell(1, 1);  % normalization (inverse of singular values)
Wm = cell(1, 1); % mask DR
aW = cell(1, 1); % preconditioner
epsilon = cell(1, 1);

H{1} = cell(realdatablocks, 1);
yT{1} = cell(realdatablocks, 1);
T{1} = cell(realdatablocks, 1);
Wm{1} = cell(realdatablocks, 1);
aW{1} = cell(realdatablocks, 1);
norm_res = cell(realdatablocks, 1);
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
    if fouRed_type == 1
        fprintf('Threshold level: remove %d percentile\n', param_fouRed.gamma);
        prob = 1 - param_fouRed.gamma / 100;
    elseif fouRed_type == 2
        fprintf('Threshold level: keep %d sigma\n', param_fouRed.gamma);
        p = normcdf([-param_fouRed.gamma param_fouRed.gamma]);
        prob = p(2) - p(1);
    end
end

precision = 1e-16;

nb_red = 0;

for j = 1:realdatablocks
    Hl{j} = (Gw{j}') * Gw{j}; % progressively write to disk? (possibly huge...)
    peak = max(max(abs(Hl{j})));
    Hl{j} = Hl{j} .* (abs(Hl{j}) > peak * precision);

    if reduction_version == 1    % Due to complex modeling and not better performance, reduction version 1 is not used any more!
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
%         Mask = (d_mat > param_fouRed.gamma);
%         th = full(peak * precision);
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

    ind_nz = d_mat > 0;       % non-zero singular values
    kept = sum(Mask(:));
    total = numel(Mask);
    percentage = kept / total * 100;
    nb_red = nb_red + kept;
    fprintf('\nBlock %d of channel %d: %d non-zero singular values, %d over %d, or %f%% of data are kept, threshold=%e\n', j, chInd, sum(ind_nz), kept, total, percentage, th);

    d_mat = d_mat(Mask);
    Hl{j} = Hl{j}(Mask, :);

    Tl{j} = d_mat;
    Tl{j} = 1 ./ sqrt(Tl{j});
    Wml{j} = Mask;

    if reduction_version == 1
        aWl{j} = 1;
        yTl{j} = dataReduce(yb{j}, Gw{j}', Wl{j}, At, Tl{j}, Wml{j});
        resRed = dataReduce(res{j}, Gw{j}', Wl{j}, At, Tl{j}, Wml{j});       % !!! for calibrated data, residuals are known
        norm_res{j} = norm(resRed);
    elseif reduction_version == 2
        aWl{j} = 1;
        Gt = Gw{j}';
        yTl{j} = Tl{j} .* (Gt(Wml{j}, :) * yb{j});
        resRed = Tl{j} .* (Gt(Wml{j}, :) * res{j});        % !!! for calibrated data, residuals are known
        norm_res{j} = norm(resRed);
        clear Gt;    % save memory
    end
    Gw{j} = [];             % save memory
    yb{j} = [];
    res{j} = [];
    fprintf('Block %d of channel %d, estimated epsilon: %f\n', j, chInd, norm_res{j});
end

fprintf('Reduced data size: %d\n', nb_red);
fprintf('H matrix memory: \n');
whos Hl;

H{1} = Hl;
W{1} = Wl;
yT{1} = yTl;
T{1} = Tl;
aW{1} = aWl;
Wm{1} = Wml;
epsilon{1} = norm_res;

DRfilename = ['/lustre/home/shared/sc004/dr_', num2str(realdatablocks), 'b_result_real_data/CYG_DR_cal_', num2str(realdatablocks), 'b_ind', ...
    num2str(subInd(1)), '_', num2str(subInd(end)), '_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), ...
    '=', num2str(chInd), '.mat'];
if ~isfile(DRfilename)
    save(DRfilename, '-v7.3', 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'epsilon', 'xsol');
end

fprintf('Dimensionality reduction and epsilon estimation are finished\n');
