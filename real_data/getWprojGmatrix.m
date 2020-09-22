function [GNUFFTW]= getWprojGmatrix(GNUFFT,wCoordinate,param,levelC,levelG)
%% input
% gImDims: image dim.
% FoV: imaged FoV in deg
% wCoordinate:  w coordinate  in UNITS OF THE WAVELENGTH (NOT RADIANS)
% [levelC,levelG]: energy thresholds for the chirp and the full kernel (should be >0.9, <1)

%% output
% GNUFFTW: G matrix incorporating the w correction
%% Global vars, Image & data dims
paddFourierFactor   = reshape(param.ox,[1,2]); % default padding factor in the Fourier domain - update otherwise
% paddFourierFactor   = param.oy;
gImDims = reshape(param.gImDims,[1,2]);
FoV     = param.FoV;
%
MeasNum = numel(wCoordinate);
uvPixelSize = 1./(paddFourierFactor.*FoV);
gGriddedFourierDim  = paddFourierFactor.*gImDims;
%% read baseline coordinates and other vars
uvHalfBW      = (uvPixelSize .*(gImDims)) ; % half imaged BW
clear UVw; % clear
eRatioWTermPrBslnx = abs(sin(FoV(1))*wCoordinate)./uvHalfBW(1); % ratio between the wterm BW and the imaged BW
eRatioWTermPrBslny = abs(sin(FoV(2))*wCoordinate)./uvHalfBW(2);
eRatioWTermPrBsln = (max(eRatioWTermPrBslnx ,eRatioWTermPrBslny));
fprintf("Eff. MAX ratio the wterm BW and the imaged BW %f\n",max(eRatioWTermPrBsln))
%% wkernel bins
% 9 bins are considered for now in the computation of the w-term,
%  more should be included for faster computations,
%  if tiny w coordinate -->  no wcorrection

eRatioWTermPrBsln = eRatioWTermPrBsln .*4;

dimBINS = 1./[64,56,48,40,32,24,16,8,4];
eRatioLowerBound = 4*max(1./max(gGriddedFourierDim));
dimBINS = dimBINS(dimBINS> eRatioLowerBound);

dimBINS =[eRatioLowerBound dimBINS];
WcoorBinPos =cell(length(dimBINS),1);
WcoorBin = zeros(length(wCoordinate),1);
WcoorBinPos{1} = find(eRatioWTermPrBsln<=eRatioLowerBound);
WcoorBin(WcoorBinPos{1})=1;
for bin = 2:length(dimBINS)-1
    WcoorBinPos{bin} = find( eRatioWTermPrBsln .*(eRatioWTermPrBsln>dimBINS(bin-1)).*(eRatioWTermPrBsln<=dimBINS(bin)));
    WcoorBin( WcoorBinPos{bin}) = bin;
end

%% w kernel
NBIS1  =  floor(min(gImDims).*dimBINS );
NBIS2  =  2*NBIS1;%floor(gImDims(2).*dimBINS );

NBIS  = [NBIS1(:)  NBIS2(:)];
nterm = cell(length(dimBINS),1);
PHASE =  cell(length(dimBINS),1);
for binW = 2 : numel(dimBINS)
    % build the lm grid
    nDim = NBIS(binW,:);
    Nt_SBIS = paddFourierFactor.* nDim;
    [l_SSBIS,m_SSBIS] = meshgrid(-Nt_SBIS(2)/2:Nt_SBIS(2)/2 -1, -Nt_SBIS(1)/2:Nt_SBIS(1)/2  -1);
    dl_SBIS = paddFourierFactor(2)*sin(FoV(1))/Nt_SBIS(1); %grid size - y-dir
    dm_SBIS = paddFourierFactor(1)*sin(FoV(2))/Nt_SBIS(2);
    l_SSBIS = l_SSBIS.*dl_SBIS;
    m_SSBIS = m_SSBIS.*dm_SBIS;
    % build the n-term (later used to build the chirp)
    nterm{binW}  = sqrt(1 - l_SSBIS.^2 - m_SSBIS.^2) ;
    nshiftSBIS =  nDim./paddFourierFactor;
    [fxSBIS,fySBIS] = meshgrid((0:(Nt_SBIS(2)-1))/Nt_SBIS(2),(0:(Nt_SBIS(1)-1))/Nt_SBIS(1).');
    omSBIS = -2*1i*pi*([fySBIS(:), fxSBIS(:)]);
    phaseSBIS = exp((omSBIS * nshiftSBIS(:))).';
    PHASE{binW}  = reshape(phaseSBIS, Nt_SBIS);
    
end
clear *BIS
%% parpool: @Ming: update numWorkers if needed
% NumWorkers = 40;
% gcp('nocreate');
% util_create_parcluster;

%% wkernel: build the w kernels
% Kernel sparsification: norm adopted
Lnorm = 2;

% prepare data for sparse conv.
gNUFFT = nnz(GNUFFT(1,:));
GNUFFT  = GNUFFT.';
[rGNUFFT1d,~,vGNUFFT] = find(GNUFFT);
clear GNUFFT; %clear
rGNUFFT1d = reshape(rGNUFFT1d,gNUFFT,MeasNum).';
vGNUFFT = reshape(vGNUFFT,gNUFFT,MeasNum).';
[rGNUFFT,cGNUFFT,~] = ind2sub(gGriddedFourierDim,rGNUFFT1d);

%init
dummyRow = cell(MeasNum,1);
dummyCol = cell(MeasNum,1);
dummyVal = cell(MeasNum,1);

% dummy vars
step = min(floor(MeasNum/2),250000);
CountVectSBIS = 1:step:MeasNum;
timeElapsed   = 0;
%
fprintf("\nINFO: Computing convolutions of the G matrix\n")
%
for row = 1:length(CountVectSBIS)-1
    wrowStart = CountVectSBIS(row);
    if row ==(length(CountVectSBIS) -1),wrowEnd = MeasNum;
    else,   wrowEnd   = CountVectSBIS(row+1) -1;
    end
    tStart  = tic;
    parfor wrow = wrowStart : wrowEnd
        %% build chirp (WKERNEL)
        WKERNELS = 1;
        binW   = WcoorBin(wrow);
        
        if  binW>1
            chirp  = exp(-2*1i*pi*wCoordinate(wrow).*(nterm{binW} -1))./nterm{binW}; % chirp analytical expr. in the image space
            chirpFourier = ifftshift(PHASE{binW}.*fft2(chirp)); % no phase is adopted in the FFT (check NUFFT shift convention)
            chirpFourier =  chirpFourier./norm(chirpFourier(:));
            % sparsification of the chirp
            [thresVal,~] =  BISECT(levelC, 1e-5, abs(chirpFourier),2);
            WKERNELS = sparse(chirpFourier.*(abs(chirpFourier)>thresVal));
            WKERNELS = WKERNELS./(norm(nonzeros(WKERNELS))); % normalisation to avoid energy loss (not compulsory)
        end
        %
        %
        %% sparse convolution[
        posOrig=[rGNUFFT(wrow,:).' cGNUFFT(wrow,:).'];
        posShifted  = shift_ind(posOrig,gGriddedFourierDim(1),gGriddedFourierDim(2));
        posOrig =[]; % clear
        kernel = [];
        kernel.i = posShifted(:,1);
        kernel.j = posShifted(:,2);
        kernel.a = vGNUFFT(wrow,:);
        kernel.dim = gGriddedFourierDim;
        posShifted = []; % clear
        
        fullKernel =  sconv2_mod_(kernel,WKERNELS,'same');
        % sparsify the G kernel
        [thresVal,~]  = BISECT(levelG, 1e-5, abs(nonzeros(fullKernel)),Lnorm);
        [dummyRowBIS,dummyColBIS,dummyVal{wrow}] = find(fullKernel.*(abs(fullKernel)>thresVal));
        
        % shift along the column dim.
        dc_Shift = shift_ind([dummyRowBIS dummyColBIS],gGriddedFourierDim(1),gGriddedFourierDim(2));
        dummyCol{wrow} = sub2ind(gGriddedFourierDim,dc_Shift(:,1),dc_Shift(:,2));
        %
        
    end
    tEnd = toc(tStart);
    fprintf("\nslice %d. Time %dmn %dsd  ",wrowStart,floor(tEnd/60),floor(mod(tEnd,60)))
    timeElapsed =timeElapsed + tEnd ;
end

clear *GNUFFT ;
%% build updated G mat
parfor wrow =1 :MeasNum
    dummyRow{wrow}= wrow*ones(size(dummyVal{wrow}));
end
tic
GNUFFTW= sparse(cell2mat(dummyRow),cell2mat(dummyCol),cell2mat(dummyVal),MeasNum,prod(gGriddedFourierDim));
toc
clear dummy*;
end


%% sconv2
function C = sconv2_mod_(A, B, shape)
% C = sconv2(A, B, shape)
% Like conv2 but suitable for convolution of sparse matrices
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date creation: 15/April/2013
% See also: conv2

if nargin < 3
    shape = 'full';
end

% [m, n] = size(A);
m = A.dim(1);
n = A.dim(2);

[p, q] = size(B);
% [i, j, a] =  find(A); %[A.i, A.j, A.a];

i = A.i ;
j = A.j ;
a = A.a ;
clear A;
[k, l, b] = find(B);

[I, K] = ndgrid(i, k);
[J, L] = ndgrid(j, l);
C = a(:)*b(:).';

switch lower(shape)
    case 'full'
        C = sparse(I(:)+K(:)-1,J(:)+L(:)-1, C(:), m+p-1, n+q-1);
    case 'valid'
        mnc = max([m-max(0,p-1),n-max(0,q-1)],0);
        i = I(:)+K(:)-p;
        j = J(:)+L(:)-q;
        b = i > 0 & i <= mnc(1) &   j > 0 & j <= mnc(2);
        C = sparse(i(b), j(b), C(b), mnc(1), mnc(2));
    case 'same'
        i = I(:)+K(:)-ceil((p+1)/2);
        j = J(:)+L(:)-ceil((q+1)/2);
        b = i > 0 & i <= m &  j > 0 & j <= n;
        C = sparse(i(b), j(b), C(b), m, n);
end

end % sconv2

%% BISECT
function[threshold,n]= BISECT(E_threshold, tolerance, abschirp,L_norm)
%calculates the threshold for keeping E_hreshold= percentile  of total
%energy of the chirp kernel, where abschirp is the absolute value of the
%row of the convoution kernel
% L_norm =1;
x0 = 0.0;
x1 = 1;
maxval = max(max(abschirp));
n = 1;
if L_norm ==2, Energy0=sqrt(sum(sum(abschirp.^2)));
else,  Energy0=(sum(sum((abschirp))));
end

Energy=1;
while((abs(x0-x1)/abs(x1))>=tolerance || (Energy-E_threshold)<0.0)
    x=x0+(x1-x0)/2.0;
    thresvalue=x*maxval;
    if L_norm ==2,Energy=sqrt(sum(sum(abschirp(abschirp>thresvalue).^2)))./Energy0;
    else,   Energy=(sum(sum((abschirp(abschirp>thresvalue)))))./Energy0;
    end
    if Energy < E_threshold,      x1=x;
    else, x0=x;
    end
    xtest=x0+(x1-x0)/2;
    thresvalue=xtest*maxval;
    if L_norm ==2,   Energy=sqrt(sum(sum(abschirp(abschirp>thresvalue).^2)))./Energy0;
    else,  Energy=(sum(sum((abschirp(abschirp>thresvalue)))))./Energy0;
    end
    n=n+1;
    if n>1000,   break;
    end
end
threshold=thresvalue;
end

function ind =shift_ind(ind,Nx,Ny)

i=ind(:,1);
j=ind(:,2);
%size assumed even
i_=zeros(length(i),1);
j_=zeros(length(j),1);
i_(i<Nx/2+1)=i(i<Nx/2+1)+Nx/2;
i_(i>Nx/2)=i(i>Nx/2)-Nx/2;
j_(j<Ny/2+1)=j(j<Ny/2+1)+Ny/2;
j_(j>Ny/2)=j(j>Ny/2)-Ny/2;

ind= [i_(:) j_(:)] ;

end