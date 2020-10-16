function func_nnls_dr_real_data_newData_wterm_part1(inputdir, outputdir, nb_aggch, chInd, subProb, reduction_version, realdatablocks, enable_klargestpercent, fouRed_gamma, fouRed_type, wterm, levelG, levelC, lowRes)

addpath ../lib/utils/
addpath ../fouRed/
addpath ../lib/operators
addpath ../lib/nufft

addpath ../../AD_build_wproj_op

speed = 299792458;  % light speed

fprintf('Channel number: %d\n', chInd);
fprintf('Subprob number: %d\n', subProb);
fprintf('Reduction version %d\n', reduction_version);
if fouRed_type == 1
    typeStr = 'perc';
elseif fouRed_type == 2
    typeStr = 'th';
end
fprintf('Data blocks: %d\n', realdatablocks);
fprintf('w projection: %d\n', wterm);
if wterm
    fprintf('Energy and chirp level: %f, %f\n', levelG, levelC);
end
fprintf('Low resolution: %d\n', lowRes);

%% Image size
if lowRes
    Nx = 2560;
    Ny = 2560;
else
    Nx = 4096;
    Ny = 4096;
end
%% Reduction
FT2 = @(x) fftshift(fft2(ifftshift(x))) / sqrt(numel(x));
% Config parameters
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
if wterm
    Kx = 7; % number of neighbours for nufft
    Ky = 7; % number of neighbours for nufft
else
    Kx = 7; % number of neighbours for nufft
    Ky = 7; % number of neighbours for nufft
end

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

if lowRes
    dl = dl/1.6;
end
fprintf("dl: %f\n", dl)

% Effective channel configuration (data-related)
% if nb_aggch == 16 || 32
bands = [33:960];
% elseif nb_aggch == 64
%     bands = [1:960];
% end

% yb = y1{1};
% nW_tmp = nW1{1};
if lowRes
    imPixelSize = 4096*1/max(Nx, Ny);    % resolution in arcsecond
else
    imPixelSize = 1;    % resolution in arcsecond
end
fprintf("Cell size: %f arc second\n", imPixelSize)
precision = 1e-16;

nb_red = 0;

% % data aggregation into effective channels
% nb_effch = floor(length(bands)/nb_aggch);
if wterm
    numworkers = 36;
    delete(gcp('nocreate'))
    cirrus_cluster = parcluster('local');   % slurm
    parpool(cirrus_cluster, numworkers, 'IdleTimeout', Inf);
end
for l = 1:length(chInd) %1:nb_effch
    
    j = chInd(l);
    
    yw = [];
    nW = [];
    uvw_wav = [];
    for k = 1:length(subProb)
        subk = subProb(k);
        if j==0
            ind = subk;
        else        
            ind = (j-1)*nb_aggch+subk+bands(1)-1;
        end
        fprintf("Effective channel: %d, channel:%d\n", j, ind)
        % uvw
        uvw1 = uvw(~final_flag{ind},:);
        uvw1(:,2) = -uvw1(:,2);   % change to coordinates [u,-v,w]
        wave = speed/Freqs(ind);
        uvw_wav = [uvw_wav; uvw1./wave];
%         % nufft
%         [~, ~, Gw, ~] = op_nufft([uvw_wav(:,2)*pi/(bmax*dl) uvw_wav(:,1)*pi/(bmax*dl)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);   % nufft using uvw in rad
        % integrate w-term into nufft
        nW = [nW; double(weights_ch(~final_flag{ind})')];
        yw = [yw; double(vis(~final_flag{ind},ind))];  
%         Gw = spdiags(nW, 0, size(uvw_wav,1), size(uvw_wav,1)) * Gw;
%         [WKERNELS,WKERNEL_NNZ]= get_wterm_wprojection(uvw_wav, Nx, Ny, imPixelSize);
    end
    fprintf("nb of points:%d\n", length(yw(:)))
    % vis
    yw = nW .* yw;    % natural-weighted vis
    
    if length(yw(:)) ~= length(nW(:)) || length(yw(:)) ~= size(uvw_wav,1) || length(nW(:)) ~= size(uvw_wav,1)
        fprintf("Error: Dimension not consistent!")
    end
    if wterm
        [G,opA,opAt,Lnorm,levelG,supports]= get_wplane_wprojection_G(uvw_wav, Nx, Ny, imPixelSize, levelG, levelC, nW, bmax, dl);
        A = opA{1};
        At = opAt{1};
    else
        [A, At, G, ~] = op_nufft([uvw_wav(:,2)*pi/(bmax*dl) uvw_wav(:,1)*pi/(bmax*dl)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);   % nufft using uvw in rad
        G = spdiags(nW, 0, size(uvw_wav,1), size(uvw_wav,1)) * G;
    end
    
    clear uvw_wav nW
    
    yw = G' * yw;   % Grided vis
    % remove zero-columns for economic RAM storage
%     Wl = G' * ones(size(G, 2), 1) ~= 0;
%     Wl = any(G, 1).';
%     G = G(:,Wl);
    
    % RAM monitoring
    fprintf('\nG with w-term matrix memory: \n')
    whos G
    
    startH = tic;
    Hl = G'*G;
    fprintf('H computation time: %e\n',toc(startH))
%         clear tmp
    fprintf('Initial H matrix memory: \n')
    whos Hl

    clear G    % save memory
    
    if wterm
        if lowRes
            Hfilename = [outputdir,'/ESO137_LOW_H_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '_subProb', num2str(subProb(1)), '_', num2str(subProb(end)), '=', num2str(j), '.mat'];
        else
            Hfilename = [outputdir,'/ESO137_H_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '_subProb', num2str(subProb(1)), '_', num2str(subProb(end)), '=', num2str(j), '.mat'];
        end
    else
        if lowRes
            Hfilename = [outputdir,'/ESO137_LOW_H_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_subProb', num2str(subProb(1)), '_', num2str(subProb(end)), '=', num2str(j), '.mat'];
        else
            Hfilename = [outputdir,'/ESO137_H_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_subProb', num2str(subProb(1)), '_', num2str(subProb(end)), '=', num2str(j), '.mat'];
        end
    end
    
    if ~isfile(Hfilename)
        save(Hfilename, '-v7.3', 'Hl', 'yw', 'A', 'At');
    end
    
%     % remove small values according to numeric precision
%     peak = max(max(abs(Hl)));
%     Hl = Hl .* (abs(Hl) > peak * precision);
%     
%     fprintf('Effective channel: %d, initial H matrix memory: \n', j)
%     whos Hl
%     
%     if reduction_version == 1    
%     %     fast matrix probing (using psf)
%         dirac2D = zeros(Ny, Nx);
%         dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;
%         PSF = operatorIpsf(dirac2D, A, At, Hl, [oy*Ny, ox*Nx]);
%         covariancemat = FT2(PSF);
%         d_mat = abs((covariancemat(:)));
%     elseif reduction_version == 2
%         d_mat = full(abs(diag(Hl)));
%     end
%     
%     if param_fouRed.enable_klargestpercent
% %         Mask = (d_mat > param_fouRed.gamma);
% %         th = full(peak * precision);
%         th = 0;
%         Mask = (d_mat > th);
% %         Mask = (d_mat >= prctile(d_mat,100-param_fouRed.klargestpercent));
%     elseif param_fouRed.enable_estimatethreshold
%     %     estimate threshold
%         d_mat_sort = sort(d_mat);
%         d_mat_sort_cumsum = cumsum(d_mat_sort);
%         d_mat_sum = d_mat_sort_cumsum(end); % whole energy of d_mat
%         th_energy = d_mat_sum * (1 - prob); % threshold according to the k-sigma rule
%         th = d_mat_sort(find(d_mat_sort_cumsum >= th_energy, 1, 'first'));
%         Mask = (d_mat >= th);
%     end
%     
%     ind_nz = d_mat > 0;       % non-zero singular values
%     kept = sum(Mask(:));
%     total = numel(Mask);
%     percentage = kept / total * 100;
%     nb_red = nb_red + kept;
%     fprintf('\nEffective channel %d: %d non-zero singular values, %d over %d, or %f%% of data are kept, threshold=%e\n', j, sum(ind_nz), kept, total, percentage, th);
%     
%     d_mat = d_mat(Mask);
%     Hl = Hl(Mask,:);
%     
%     fprintf('Reduced H matrix memory: \n')
%     whos Hl
% 
%     Tl = d_mat;
%     Tl = 1./sqrt(Tl);
%     Wml = Mask;
%     
%     clear d_mat Mask
% %     H{1}{1} = Hl;
% %     W{1}{1} = Wl;
% %     T{1}{1} = Tl;
% %     Wm{1}{1} = Wml;
% 
% %     Hfilename = [outputdir,'/ESO137_H_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '=', num2str(chInd), '.mat'];
% %     if ~isfile(Hfilename)
% %         save(Hfilename, '-v7.3', 'H', 'W', 'T', 'Wm');
% %     end
%     
%     if reduction_version == 1
%         aWl = 1;
%         im = FT2(real(At(yw)));
%         im = im(:);
%         yTl = Tl .* im(Wml);
% %         resRed = dataReduce(res{j}, Gw{j}', Wl{j}, At, Tl{j}, Wml{j});       % !!! for calibrated data, residuals are known
% %         norm_res{j} = norm(resRed);
%     elseif reduction_version == 2
%         aWl = 1;
%         yw = yw(Wl,:);
%         yTl = Tl.*yw(Wml,:);
% %         resRed = Tl{j}.*(Gt(Wml{j},:) * res{j});        % !!! for calibrated data, residuals are known
% %         norm_res{j} = norm(resRed);
%         eps_normy = 0.01 * norm(yTl(:));
%         fprintf("Epsilon estimated from norm(y): %e\n", eps_normy)
%         [~, norm_res] = fb_dr_nnls(yTl, A, At, Hl, Wl, Tl, Wml, param_nnls, reduction_version);
%     end
%                  % save memory
%     fprintf('Effective channel %d, estimated epsilon: %f\n', j, norm_res)
% 
%     H{1}{1} = Hl;
%     W{1}{1} = Wl;
%     yT{1}{1} = yTl;
%     T{1}{1} = Tl;
%     aW{1}{1} = aWl;
%     Wm{1}{1} = Wml;
%     epsilon{1}{1} = norm_res;
%     epsilon_normy{1}{1} = eps_normy;
% 
%     DRfilename = [outputdir,'/ESO137_DR_fouRed', num2str(reduction_version), '_', typeStr, num2str(fouRed_gamma), '_w', num2str(levelG), '_', num2str(levelC), '=', num2str(chInd), '.mat'];
%     if ~isfile(DRfilename)
%         save(DRfilename, '-v7.3', 'A','At','H', 'W', 'T', 'Wm', 'yT', 'aW', 'epsilon', 'epsilon_normy');
%     end

end

% fprintf('Reduced data size: %d\n', nb_red)

fprintf('Dimensionality reduction and epsilon estimation part 1/2 are finished\n')

