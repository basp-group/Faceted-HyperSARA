function func_sum_H_dr(inputdir, outputdir, chInd, nbAgg, group, reduction_version, realdatablocks, enable_klargestpercent, fouRed_gamma, fouRed_type)

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

Nx = 4096;
Ny = 4096;

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

% parameter NNLS
param_nnls.verbose = 2;       % print log or not
param_nnls.rel_obj = 1e-4;    % stopping criterion
param_nnls.max_iter = 2000;     % max number of iterations 1000
param_nnls.sol_steps = [inf]; % saves images at the given iterations
param_nnls.beta = 1;

fprintf('Dimensionality reduction begins ...\n')
fprintf('Enable k-largest percentage: %d\n', param_fouRed.enable_klargestpercent);
fprintf('Enable automatic threshold: %d\n', param_fouRed.enable_estimatethreshold);

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

[A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);

Hl1 = sparse(ox*Nx*oy*Ny,ox*Nx*oy*Ny);
yl1 = zeros(ox*Nx*oy*Ny,1);
for i = 1:group:nbAgg
    filename = [inputdir, 'ESO137_H_fouRed2_perc5_subProb', num2str(i), '_', num2str(i+group-1), '=', num2str(chInd), '.mat'];
    fprintf('Read file name: %s\n', filename)
    load(filename, 'Hl', 'yw')
    Hl1 = Hl1 + Hl;
    clear Hl
    yl1 = yl1 + yw;
    clear yw
end
yw = yl1;
% remove zero-columns for economic RAM storage
Wl = Hl1 * ones(size(Hl1, 1), 1) ~= 0;
Hl = Hl1(Wl, Wl);
clear Hl1

precision = 1e-16;
% remove small values according to numeric precision
peak = max(max(abs(Hl)));
Hl = Hl .* (abs(Hl) > peak * precision);

fprintf('Effective channel: %d, initial H matrix memory: \n', chInd)
whos Hl

nb_red = 0;
if reduction_version == 1    
%     fast matrix probing (using psf)
    dirac2D = zeros(Ny, Nx);
    dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;
    PSF = operatorIpsf(dirac2D, A, At, Hl, [oy*Ny, ox*Nx]);
    covariancemat = FT2(PSF);
    d_mat = abs((covariancemat(:)));
elseif reduction_version == 2
    d_mat = full(abs(diag(Hl)));
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
fprintf('\nEffective channel %d: %d non-zero singular values, %d over %d, or %f%% of data are kept, threshold=%e\n', chInd, sum(ind_nz), kept, total, percentage, th);

d_mat = d_mat(Mask);
Hl = Hl(Mask,:);

fprintf('Reduced H matrix memory: \n')
whos Hl

Tl = d_mat;
Tl = 1./sqrt(Tl);
Wml = Mask;

% H{1}{1} = Hl;
% W{1}{1} = Wl;
% T{1}{1} = Tl;
% Wm{1}{1} = Wml;

% Hfilename = [outputdir,'/ESO137_H_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(chInd), '.mat'];
% if ~isfile(Hfilename)
%     save(Hfilename, '-v7.3', 'H', 'W', 'T', 'Wm');
% end

if reduction_version == 1
    aWl = 1;
    im = FT2(real(At(yw)));
    im = im(:);
    yTl = Tl .* im(Wml);
%         resRed = dataReduce(res{j}, Gw{j}', Wl{j}, At, Tl{j}, Wml{j});       % !!! for calibrated data, residuals are known
%         norm_res{j} = norm(resRed);
elseif reduction_version == 2
    aWl = 1;
    yw = yw(Wl,:);
    yTl = Tl.*yw(Wml,:);
%         resRed = Tl{j}.*(Gt(Wml{j},:) * res{j});        % !!! for calibrated data, residuals are known
%         norm_res{j} = norm(resRed);
    [~, norm_res] = fb_dr_nnls(yTl, A, At, Hl, Wl, Tl, Wml, param_nnls, reduction_version);
end
             % save memory
fprintf('Effective channel %d, estimated epsilon: %f\n', chInd, norm_res)

H{1}{1} = Hl;
W{1}{1} = Wl;
yT{1}{1} = yTl;
T{1}{1} = Tl;
aW{1}{1} = aWl;
Wm{1}{1} = Wml;
epsilon{1}{1} = norm_res;

DRfilename = [outputdir,'/ESO137_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(chInd), '.mat'];
if ~isfile(DRfilename)
    save(DRfilename, '-v7.3', 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'epsilon');
end

fprintf('Reduced data size: %d\n', nb_red)

fprintf('Dimensionality reduction and epsilon estimation are finished\n')


