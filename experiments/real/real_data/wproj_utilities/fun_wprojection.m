function Gw = fun_wprojection(GNUFFT, W, flag, param)
%% input
% gImDims: image dim.
% FoV: imaged FoV in deg
% wCoordinate:  w coordinate  in UNITS OF THE WAVELENGTH (NOT RADIANS)
% [levelC,levelG]: energy thresholds for the chirp and the full kernel (should be >0.9, <1)

%% output
% GNUFFTW: G matrix incorporating the w correction
% supports: kernel dimensions

%% Global vars, Image & data dims
paddFourierFactor   = param.paddFourierFactor; % default padding factor in the Fourier domain - update otherwise
CEnergyL2 = param.CEnergyL2;
GEnergyL2 = param.GEnergyL2;
gImDims = param.gImDims;
FoV     = [param.FoVy param.FoVx];
% ProjBaseline = param.ProjBaseline; %(sqrt(u.^2+v.^2))
MeasNum = numel(W);
uvPixelSize = 1 ./ (paddFourierFactor .* FoV);
ImFourierDims  = paddFourierFactor .* gImDims;
% read baseline coordinates and other vars
uvHalfBW      = [param.vGridSize param.uGridSize] .* gImDims; % half imaged BW

%% wkernel bins
[WResolutionBin, PHASE, nTerm] = get_chirp_details(W, FoV, gImDims, paddFourierFactor, uvHalfBW);

%% wkernel: build the w kernels
% Kernel sparsification: norm adopted
Lnorm = 2;

% prepare data for sparse conv.
if ~isempty(flag); nActiveRows = nnz(flag);
else;  nActiveRows = numel(W);
end
gNUFFT = param.supportK;
% GNUFFT  = GNUFFT.';
[rGNUFFT1d, ~, vGNUFFT] = find(GNUFFT);
clear GNUFFT; % clear
rGNUFFT1d = reshape(rGNUFFT1d, gNUFFT, nActiveRows).';
vGNUFFT   = reshape(vGNUFFT, gNUFFT, nActiveRows).';
[rGNUFFT, cGNUFFT, ~] = ind2sub(ImFourierDims, rGNUFFT1d);
rGNUFFT1d = [];
rGNUFFT = rGNUFFT.';
cGNUFFT = cGNUFFT.';
vGNUFFT = vGNUFFT.';
% init
tStart  = tic;
if ~isempty(flag)
    wCoordinate_flagged  = W(flag > 0); clear wCoordinate;
    WcoorBin_flagged =  WResolutionBin(flag > 0); clear WcoorBin;
    ActiveRows = find(flag);
else
    wCoordinate_flagged = W; clear wCoordinate;
    WcoorBin_flagged = WResolutionBin; clear WcoorBin;
    ActiveRows = ':';
end
dummyCol = cell(nActiveRows, 1);
dummyVal = cell(nActiveRows, 1);
dummyRow = cell(nActiveRows, 1);
SupportGKernel = zeros(nActiveRows, 1);
parfor wrow = 1:nActiveRows
    %% build chirp (WKERNEL)
    w_kernel = 1;
    if  WcoorBin_flagged(wrow) > 1
        w_kernel = get_chirp(W(wrow), nTerm{WcoorBin_flagged(wrow)}, PHASE{WcoorBin_flagged(wrow)}, CEnergyL2);

        %% sparse convolution[
        nufft_kernel = [];
        nufft_kernel.dim = ImFourierDims;
        nufft_kernel.i = rGNUFFT(:, wrow);
        nufft_kernel.j = cGNUFFT(:, wrow);
        nufft_kernel.a =  vGNUFFT(:, wrow);
        % conv
        full_kernel = (sconv2_modified(nufft_kernel, w_kernel, 'same')); % conj for G^\daggerger
    else
        full_kernel = sparse(rGNUFFT(:, wrow), cGNUFFT(:, wrow), vGNUFFT(:, wrow), ImFourierDims(1), ImFourierDims(2));
    end
    % sparsify the G kernel

    if GEnergyL2 < 1
        [thresVal, ~] = BISECT(GEnergyL2, 1e-5, abs(nonzeros(full_kernel)), Lnorm);
        full_kernel = full_kernel .* (abs(full_kernel) > thresVal);
    end
    SupportGKernel(wrow) = nnz(full_kernel);
    % ifftshift
    [posShift1, posShift2, dummyVal{wrow}] = find(conj(full_kernel));
    posOrig = shift_ind([posShift1, posShift2], ImFourierDims(1), ImFourierDims(2));
    %
    dummyCol{wrow} = sub2ind(ImFourierDims, posOrig(:, 1), posOrig(:, 2));

    if strcmp(ActiveRows, ':');  activeRowId = wrow;
    else; activeRowId = ActiveRows(wrow);
    end

    dummyRow{wrow} = activeRowId * ones(SupportGKernel(wrow), 1);

end
timeElapsed = toc(tStart);
% fprintf('\nkernels computed %dmn \n ',floor(timeElapsed/60));
clear *GNUFFT nufft_kernel full_kernel w_kernel;
%% build updated G mat

if ~isempty(flag)
    KernelG_Support = zeros(MeasNum, 1);
    KernelG_Support(flag > 0) = SupportGKernel;
else; KernelG_Support = SupportGKernel;
end
clear SupportGKernel;
dummyCol = cell2mat(dummyCol);
dummyRow = cell2mat(dummyRow);
dummyVal = cell2mat(dummyVal);

% fprintf('\nBuilding conj(transpose(G)) for mem reasons  ... ');
FourierDim = prod(ImFourierDims);

Gw = sparse(dummyCol, dummyRow, dummyVal, FourierDim, MeasNum);
clear dummy*;
KernelG_activeMode = []; % find(sum(abs(GNUFFTW),2));

end

% eRatioWTermPrBsln = abs(sin(FoV)*wCoordinate)./uvHalfBW; % ratio between the wterm BW and the imaged BW
% eResolution = uvHalfBW/max(ProjBaseline); % ratio between the imaged BW and the probed BW
% % fprintf("\nEff. ratio  the imaged BW and the probed BW %f\n",eResolution )
% % fprintf("Eff. MAX ratio the wterm BW and the imaged BW %f\n",max(eRatioWTermPrBsln))
