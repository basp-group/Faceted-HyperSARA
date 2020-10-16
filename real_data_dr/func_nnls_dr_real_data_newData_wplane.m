function func_nnls_dr_real_data_newData_wplane(inputdir, outputdir, nb_aggch, chInd, reduction_version, realdatablocks, enable_klargestpercent, fouRed_gamma, fouRed_type)

addpath ../lib/utils/
addpath ../fouRed/
addpath ../lib/operators
addpath ../lib/nufft

addpath ../../AD_build_wproj_wplanes_op

speed = 299792458;  % light speed

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

% instantiate variables for DR
% H = cell(1, 1);  % holographic matrices
% yT = cell(1, 1); % reduced data
% T = cell(1, 1);  % normalization (inverse of singular values)
% Wm = cell(1, 1); % mask DR
% aW = cell(1, 1); % preconditioner
% epsilon = cell(1, 1);
% 
% H{1} = cell(realdatablocks, 1);
% yT{1} = cell(realdatablocks, 1);
% T{1} = cell(realdatablocks, 1);
% Wm{1} = cell(realdatablocks, 1);
% aW{1} = cell(realdatablocks, 1);
% norm_res = cell(realdatablocks, 1);
% % sing{1} = cell(numel(y_tmp), 1);
% 
% Hl = H{1};
% yTl = yT{1};
% Tl = T{1};
% Wml = Wm{1};
% aWl = aW{1};

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

% load data
datafilename = [inputdir,'ESO137_FULL.mat'];
fprintf("Read file: %s\n", datafilename)
load(datafilename,'weights_ch','vis','Freqs','uvw');
% load final flag
flagfilename = [inputdir,'final_flag.mat'];
fprintf("Read file: %s\n", flagfilename)
load(flagfilename, 'bmax', 'dl', 'final_flag')

% Effective channel configuration (data-related)
if nb_aggch == 16
    bands = [33:960];
elseif nb_aggch == 64
    bands = [1:960];
end

% yb = y1{1};
% nW_tmp = nW1{1};
imPixelSize = 1;    % resolution in arcsecond
precision = 1e-16;

nb_red = 0;

% % data aggregation into effective channels
% nb_effch = floor(length(bands)/nb_aggch);
numworkers = 36;
for l = 1:length(chInd) %1:nb_effch
    
    j = chInd(l);
    
    Hl = sparse(ox*Nx*oy*Ny,ox*Nx*oy*Ny);   % Holographic matrix
    yw = 0;
    for k = 1:nb_aggch
        ind = (j-1)*nb_aggch+k+bands(1)-1;
        fprintf("Effective channel: %d, channel:%d\n", j, ind)
        % uvw
        uvw_wav = uvw(~final_flag{ind},:);
        uvw_wav(:,2) = -uvw_wav(:,2);   % change to coordinates [u,-v,w]
        wave = speed/Freqs(ind);
        uvw_wav = uvw_wav./wave;
        % nufft
        [~, ~, Gw, ~] = op_nufft([uvw_wav(:,2)*pi/(bmax*dl) uvw_wav(:,1)*pi/(bmax*dl)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);   % nufft using uvw in rad
        % integrate w-term into nufft
        nW = double(weights_ch(~final_flag{ind})');
        Gw = spdiags(nW, 0, size(uvw_wav,1), size(uvw_wav,1)) * Gw;
        [WKERNELS,WKERNEL_NNZ,coefA,BlocksIdx]= get_wplane_wprojection(uvw_wav,Nx,Ny,imPixelSize);
        G = update_G_wprojection(Gw, WKERNELS, WKERNEL_NNZ,[Ny Nx]);
        % vis
        y = double(vis(~final_flag{ind},ind));
        nW = double(weights_ch(~final_flag{ind})');
        yw = yw + Gw'*(nW.*y);      % Gridded natural-weighted vis
        % RAM monitoring
        fprintf('G matrix memory: \n')
        whos G
        % split G into smaller Gs to accelerate
        Htmp = sparse(ox*Nx*oy*Ny,ox*Nx*oy*Ny);
        nbRows = size(G,1);
        slices = round(nbRows / numworkers);
        parfor m = 1:numworkers
            ind1 = (m-1)*slices + 1;
            ind2 = m * slices;
            if ind2 > nbRows
                ind2 = nbRows;
            end
            Htmp = Htmp + Gw(ind1:ind2,:)'* G(ind1:ind2,:);             % H = G0' * G;
        end
        fprintf('Initial H matrix memory: \n')
        whos Htmp
        % Add splitted H together
        Hl = Hl + Htmp;
        clear uvw_wav y nW Gw G Htmp    % save memory
    end
%     % remove zero-columns for economic RAM storage
%     Wl = Hl * ones(size(Hl, 1), 1) ~= 0;
%     Hl = Hl(Wl, Wl);
    % remove small values according to numeric precision
    peak = max(max(abs(Hl)));
    Hl = Hl .* (abs(Hl) > peak * precision);
    
    fprintf('Effective channel: %d, initial H matrix memory: \n', j)
    whos Hl
    
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
    fprintf('\nEffective channel %d: %d non-zero singular values, %d over %d, or %f%% of data are kept, threshold=%e\n', j, sum(ind_nz), kept, total, percentage, th);
    
    d_mat = d_mat(Mask);
    Hl = Hl(Mask,:);
    
    fprintf('Reduced H matrix memory: \n')
    whos Hl

    Tl = d_mat;
    Tl = 1./sqrt(Tl);
    Wml = Mask;
    
    if reduction_version == 1
        aWl = 1;
        im = FT2(real(At(yw)));
        im = im(:);
        yTl = Tl .* im(Wml);
%         resRed = dataReduce(res{j}, Gw{j}', Wl{j}, At, Tl{j}, Wml{j});       % !!! for calibrated data, residuals are known
%         norm_res{j} = norm(resRed);
    elseif reduction_version == 2
        aWl = 1;
        yTl = Tl.*yw(Wml,:);
%         resRed = Tl{j}.*(Gt(Wml{j},:) * res{j});        % !!! for calibrated data, residuals are known
%         norm_res{j} = norm(resRed);
        [~, norm_res] = fb_dr_nnls(yTl, A, At, Hl, Wl, Tl, Wml, param_nnls, reduction_version);
    end
                 % save memory
    fprintf('Effective channel %d, estimated epsilon: %f\n', j, norm_res)

    H{1}{1} = Hl;
    W{1}{1} = [];
    yT{1}{1} = yTl;
    T{1}{1} = Tl;
    aW{1}{1} = aWl;
    Wm{1}{1} = Wml;
    epsilon{1}{1} = norm_res;

    DRfilename = [outputdir,'/ESO137_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(chInd), '.mat'];
    if ~isfile(DRfilename)
        save(DRfilename, '-v7.3', 'H', 'W', 'yT', 'T', 'aW', 'Wm', 'epsilon');
    end

end

fprintf('Reduced data size: %d\n', nb_red)

fprintf('Dimensionality reduction and epsilon estimation are finished\n')

