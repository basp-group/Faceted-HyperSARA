function func_nnls_dr_real_data_newData(inputdir,outputdir,chInd, reduction_version, realdatablocks, enable_klargestpercent, fouRed_gamma, fouRed_type)

addpath ../lib/utils/
addpath ../fouRed/
addpath ../lib/operators
addpath ../lib/nufft

fprintf('Channel number: %d\n', chInd);
fprintf('Reduction version %d\n', reduction_version);
if fouRed_type == 1
    typeStr = 'perc';
elseif fouRed_type == 2
    typeStr = 'th';
end
fprintf('Data blocks: %d\n', realdatablocks);

%% Image size
Nx = 4096;
Ny = 4096;

%% Reduction
FT2 = @(x) fftshift(fft2(ifftshift(x))) / sqrt(numel(x));
% Config parameters
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

% Fourier reduction parameters
param_fouRed.enable_klargestpercent = enable_klargestpercent;
param_fouRed.enable_estimatethreshold = ~enable_klargestpercent;
param_fouRed.gamma = fouRed_gamma;  % reduction parameter
param_fouRed.diagthresholdepsilon = 0;
param_fouRed.fastCov = 1;

% parameter NNLS
param_nnls.verbose = 2;       % print log or not
param_nnls.rel_obj = 1e-4;    % stopping criterion
param_nnls.max_iter = 2000;     % max number of iterations 1000
param_nnls.sol_steps = [inf]; % saves images at the given iterations
param_nnls.beta = 1;

fprintf('Dimensionality reduction begins ...\n')
fprintf('Enable k-largest percentage: %d\n', param_fouRed.enable_klargestpercent);
fprintf('Enable automatic threshold: %d\n', param_fouRed.enable_estimatethreshold);

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
        prob = 1 - param_fouRed.gamma/100;
    elseif fouRed_type == 2
        fprintf('Threshold level: keep %d sigma\n', param_fouRed.gamma);
        p = normcdf([-param_fouRed.gamma param_fouRed.gamma]);
        prob = p(2) - p(1);
    end
end

datafilename = [inputdir,'/data_ESO137_full_ch', num2str(chInd), '.mat'];
fprintf("Read file: %s\n", datafilename)
load(datafilename);

u_tmp = u1{1};
v_tmp = v1{1};
yb = y1{1};
nW_tmp = nW1{1};
natWy = true;

precision = 1e-16;

nb_red = 0;

for j = 1:realdatablocks
    
    [A, At, Gw, Wl] = op_p_nufft([v_tmp(j) u_tmp(j)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW_tmp(j));
    
%     Wl{1}= Gw{1}' * ones(size(Gw{1}, 1), 1) ~= 0;
%     Gw{1} = Gw{1}(:, Wl{1});
    
    fprintf('G matrix memory: \n')
    whos Gw
    
    Hl{j} = (Gw{1}')*Gw{1}; % progressively write to disk? (possibly huge...)
    fprintf('Initial H matrix memory: \n')
    whos Hl
    
    peak = max(max(abs(Hl{j})));
    Hl{j} = Hl{j} .* (abs(Hl{j}) > peak * precision);
    
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
    Hl{j} = Hl{j}(Mask,:);

    Tl{j} = d_mat;
    Tl{j} = 1./sqrt(Tl{j});
    Wml{j} = Mask;
    
%     if natWy
        fprintf('Data natural-weighted\n')
        yb{j} = yb{j} .* nW_tmp{j};
%     end
    
    if reduction_version == 1
        aWl{j} = 1;
        yTl{j} = dataReduce(yb{j}, Gw{j}', Wl{j}, At, Tl{j}, Wml{j});
%         resRed = dataReduce(res{j}, Gw{j}', Wl{j}, At, Tl{j}, Wml{j});       % !!! for calibrated data, residuals are known
%         norm_res{j} = norm(resRed);
    elseif reduction_version == 2
        aWl{j} = 1;
        Gt = Gw{j}';
        yTl{j} = Tl{j}.*(Gt(Wml{j},:) * yb{j});
%         resRed = Tl{j}.*(Gt(Wml{j},:) * res{j});        % !!! for calibrated data, residuals are known
%         norm_res{j} = norm(resRed);
        norm_res{j} = 0.01*norm(yTl{j});
%         [~, norm_res{j}] = fb_dr_nnls(yTl{j}, A, At, Hl{j}, Wl{j}, Tl{j}, Wml{j}, param_nnls, reduction_version);
        clear Gt    % save memory
    end
    Gw{j} = [];             % save memory
    yb{j} = [];
    fprintf('Block %d of channel %d, estimated epsilon: %f\n', j, chInd, norm_res{j})
end

fprintf('Reduced data size: %d\n', nb_red)
fprintf('H matrix memory: \n')
whos Hl

H{1} = Hl;
W{1} = Wl;
yT{1} = yTl;
T{1} = Tl;
aW{1} = aWl;
Wm{1} = Wml;
epsilon{1} = norm_res;

DRfilename = [outputdir,'/ESO137_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(chInd), '.mat'];
if ~isfile(DRfilename)
    save(DRfilename, '-v7.3', 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'epsilon');
end

fprintf('Dimensionality reduction and epsilon estimation are finished\n')

