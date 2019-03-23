function [Ipsf, Mask, d12, FIpsf, FIpsf_t] = fourierReduction(Gw, A, At, imsize, W, param)

% Flags monitoring

if ~isfield(param,'enable_klargestpercent') && ~isfield(param,'enable_estimatethreshold')
    param.enable_klargestpercent = 1;
end

if param.enable_klargestpercent
    if ~isfield(param, 'klargestpercent') 
        param.klargestpercent = 100; 
    end
    klargestpercent = param.klargestpercent;
elseif param.enable_estimatethreshold
    if ~isfield(param, 'gamma') 
        param.gamma = 3; 
    end
    if ~isfield(param, 'x2') && ~isfield(param, 'sigma_noise')
        error('Either ||x||_2 or sigma_noise is missing for the estimation of the threshold');
    end
    sigma_noise = param.sigma_noise;
    gamma = param.gamma;
    x2 = param.x2;
    dirty2 = param.dirty2;
end

% Flag to pull up the values of elements of the holographic matrix
% This is to avoid having VERY small values which might later explode
% during computation of inverse or reciprocal.
if ~isfield(param, 'diagthresholdepsilon') 
    param.diagthresholdepsilon = 1e-10; 
end

if ~isfield(param, 'covmatfileexists')
    param.covmatfileexists = 0;
end

if ~isfield(param, 'covmatfile')
    param.covmatfile = 'covariancemat.mat';
end

% Fast covariance matrix computation
if ~isfield(param, 'fastCov')
    param.fastCov = 1;
end

diagthresholdepsilon = param.diagthresholdepsilon;
covmatfileexists = param.covmatfileexists;
covmatfile = param.covmatfile;
fastCov = param.fastCov;

Ny = imsize(1);
Nx = imsize(2);
% Compute holographic matrix
H = Gw'*Gw;
% Create the PSF operator
Ipsf = @(x) operatorIpsf(x, A, At, H, [Ny, Nx], W);     % Phi^T Phi = At G' G A = At H A: image -> image

serialise = @(x) x(:);

fprintf('\nComputing covariance matrix...');
% Covariance operator F Phi^T Phi F^T = F At H A F^T = F Ipsf F^T
FT2 = @(x) fftshift(fft2(ifftshift(x)));
IFT2 = @(x) fftshift(ifft2(ifftshift(x)));

covoperator = @(x) serialise(FT2(Ipsf(real(IFT2(reshape(x, Ny, Nx))))));
diagonly = 1; % Only compute diagonal of covariance matrix FPhi^TPhiF^T
if covmatfileexists
    fprintf('\nLoading covariance matrix from file...');
    load(covmatfile, 'covariancemat');
else
    tstartcovmat = tic;
    if fastCov
        dirac2D = zeros(Ny, Nx);
        dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;

        PSF = Ipsf(dirac2D);
        covariancemat = FT2(PSF);
    else
        covariancemat = guessmatrix(diagonly, covoperator, Ny*Nx, Ny*Nx);
    end
    fprintf('\nSaving covariance matrix...\n');
    save(covmatfile, 'covariancemat');
    tendcovmat = toc(tstartcovmat);
    fprintf('Time to compute covariance matrix: %e s\n', tendcovmat)
end

if fastCov
%     d = abs(covariancemat(:));
    d = abs(real(covariancemat(:)));
else
    d = diag(covariancemat); %*(sigma_noise^2)
%     d = abs(d);
    d = abs(real(d));
end

% Singular values thresholding
fprintf('\nPruning covariancemat according to eigenvalues (diagonal elements)...\n');
if param.enable_klargestpercent
    Mask = (d >= prctile(d,100-klargestpercent));
elseif param.enable_estimatethreshold
    % Embed the noise
    noise = sigma_noise/sqrt(2) * (randn(size(Gw, 1),1) + 1j * randn(size(Gw, 1), 1));
    rn = FT2(At(Gw'*noise));  % Apply F Phi
    th = gamma * std(rn(:)) / x2;
%     th_dirty = gamma * std(rn(:)) / dirty2;
%     param.maxd = max(d);
%     param.th = th;
%     param.th_dirty = th_dirty;
%     param.singInit = d;
    fprintf('\nThe estimate threshold using ground truth is %e \n', th);
%     fprintf('\nThe estimate threshold using dirty image is %e \n', th_dirty);
    Mask = (d >= th_dirty);
end
d = d(Mask);
fprintf('\nThe threshold is %e \n', min(d));

d = max(diagthresholdepsilon, d);  % This ensures that inverting the values will not explode in computation
d12 = 1./sqrt(d);

FIpsf = @(x) serialise(FT2(Ipsf(x)));  % F * Ipsf, image -> vect
FIpsf_t = @(x) Ipsf(IFT2(reshape(full(x), imsize)));  % Ipsf * F^T, vect -> image

end

function y = operatorIpsf(x, A, At, H, No, W)
Ny = No(1);
Nx = No(2);

x1 = A(real(x));
if exist('W', 'var') && exist('No', 'var')
    x1 = x1(W);
end
x2 = H * x1;
if exist('W', 'var') && exist('No', 'var')
    x3 = zeros(Ny * Nx, 1);
    x3(W) = x2;
    x2 = x3;
end
y = real(At(x2));
end

    
