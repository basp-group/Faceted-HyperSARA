function [R, Rt, d12, Mask] = fourierReduction(Gw, A, At, imsize, param)
% Flags

% Flag to pull up the values of elements of the holographic matrix
% This is to avoid having VERY small values which might later explode
% during computation of inverse or reciprocal.
if ~isfield(param, 'klargestpercent') 
    param.klargestpercent = 100; 
end

if ~isfield(param, 'diagthresholdepsilon') 
    param.diagthresholdepsilon = 1e-10; 
end

if ~isfield(param, 'covmatfileexists')
    param.covmatfileexists = 0;
end

if ~isfield(param, 'covmatfile')
    param.covmatfile = 'covariancemat.mat';
end

klargestpercent = param.klargestpercent;
diagthresholdepsilon = param.diagthresholdepsilon;
covmatfileexists = param.covmatfileexists;
covmatfile = param.covmatfile;


% Flag to load from a previously saved covariance matrix file
% covmatfileexists = 0;       % Read precomputed matrix 
% covmatfile = 'data/savedfiles/covariancemat.mat';

% Flag to set if we want to approximate D with an
% identity matrix. Reset the flag to use the normal D.
% (D is the diagonal aproximation of the covariance matrix)
Ny = imsize(1);
Nx = imsize(2);
%     % Compute holographic matrix
%     h = Gw'*Gw;

% Create the new measurement operator
serialise = @(x) x(:);
R = @(x) serialise(fft2(real(At(Gw'*x(:)))));   % F Phi^T: vect -> vect
Rt = @(x) Gw*serialise(A(real(ifft2(reshape(full(x), Ny, Nx)))));   % Phi F^T: vect -> vect

%     grid2img_fwd = @(x) At(h*A(x)); % Phi^TPhi; input = [Ny, Nx] image; output = [Ny, Nx] matrix
%     grid2img_adj = @(x) At(h*A(reshape(x, Ny, Nx))); % input = Ny*Nx vector; output = [Ny, Nx] matrix
fprintf('\nComputing covariance matrix...');
% Takes a vectorized input
%     covoperator = @(x) serialise(fft2(grid2img_fwd((Ny*Nx)*ifft2(reshape(full(x), [Ny, Nx])))));
covoperator = @(x) R(Rt(x));
diagonly = 1; % Only compute diagonal of covariance matrix FPhi^TPhiF^T
if covmatfileexists
    fprintf('\nLoading covariance matrix from file...');
    load(covmatfile, 'covariancemat');
else
    tstartcovmat = tic;
    covariancemat = guessmatrix(diagonly, covoperator, Ny*Nx, Ny*Nx);
    fprintf('\nSaving covariance matrix...\n');
    save(covmatfile, 'covariancemat');
    tendcovmat = toc(tstartcovmat);
    fprintf('Time to compute covariance matrix: %e s\n', tendcovmat)
end

d = diag(covariancemat); %*(sigma_noise^2)
d = abs(d);

% d = ones(size(d)); % Disable weighting, simply do FPhi^TPhiF^T. Throw away the diagonal of covariancemat
fprintf('\nPruning covariancemat according to eigenvalues (diagonal elements)...\n');
Mask = (d >= prctile(d,100-klargestpercent));
d = d(Mask);
d = max(diagthresholdepsilon, d);  % This ensures that inverting the values will not explode in computation
d12 = 1./sqrt(d);
% Mask = sparse(1:length(nonzerocols), nonzerocols, ones(length(nonzerocols), 1), length(nonzerocols), (Ny*Nx));

