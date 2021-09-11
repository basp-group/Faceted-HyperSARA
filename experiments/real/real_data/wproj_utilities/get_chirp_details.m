function [WcoorBin, PHASE, nterm] = get_chirp_details(W, FoV, gImDims, paddFourierFactor, uvHalfBW)

%  if tiny w coordinate -->  no wcorrection
ImFourierDim = paddFourierFactor .* gImDims;
eRatioWTermPrBsln = max(abs(sin(FoV) .* W) ./ uvHalfBW); % ratio between the wterm BW and the imaged BW
% % % % eResolution = uvHalfBW/max(ProjBaseline); % ratio between the imaged BW and the probed BW
% % % % fprintf("\nEff. ratio  the imaged BW and the probed BW %f\n",eResolution )
% fprintf("Eff. MAX ratio the wterm BW and the imaged BW %f\n",max(eRatioWTermPrBsln))

eRatioWTermPrBsln = eRatioWTermPrBsln .* 5;
dimBINS = 1 ./ [256:-8:8 4];
eRatioLowerBound = 4 * max(1 ./ max(ImFourierDim));
dimBINS = dimBINS(dimBINS > eRatioLowerBound);

dimBINS = [eRatioLowerBound dimBINS];
WcoorBinPos = cell(length(dimBINS), 1);
WcoorBin = zeros(length(W), 1);
WcoorBinPos{1} = find(eRatioWTermPrBsln <= eRatioLowerBound);
WcoorBin(WcoorBinPos{1}) = 1;
for bin = 2:length(dimBINS) - 1
    WcoorBinPos{bin} = find(eRatioWTermPrBsln .* (eRatioWTermPrBsln > dimBINS(bin - 1)) .* (eRatioWTermPrBsln <= dimBINS(bin)));
    WcoorBin(WcoorBinPos{bin}) = bin;
end

%% w kernel in the image domain
NCurr_1  =  (floor(gImDims(1) .* dimBINS) - mod(floor(gImDims(1) .* dimBINS), 2));
NCurr_2  =  (floor(gImDims(2) .* dimBINS) - mod(floor(gImDims(2) .* dimBINS), 2));

NCurr  = [NCurr_1(:)  NCurr_2(:)];
nterm = cell(length(dimBINS), 1);
PHASE =  cell(length(dimBINS), 1);
for binW = 1:numel(dimBINS)
    % build the lm grid
    Nt_SCurr =  paddFourierFactor .* NCurr(binW, :);
    [l_SSCurr, m_SSCurr] = meshgrid(-Nt_SCurr(1) / 2:Nt_SCurr(1) / 2 - 1, -Nt_SCurr(1) / 2:Nt_SCurr(2) / 2  - 1);
    dl_SCurr = paddFourierFactor(1) .* sin(FoV(1)) / Nt_SCurr(1); % grid size - y-dir
    dm_SCurr = paddFourierFactor(2) .* sin(FoV(2)) / Nt_SCurr(2);
    l_SSCurr = l_SSCurr .* dl_SCurr;
    m_SSCurr = m_SSCurr .* dm_SCurr;
    % build the n-term (later used to build the chirp)
    nterm{binW}  = sqrt(1 - l_SSCurr.^2 - m_SSCurr.^2) - 1;
    nshiftSCurr =  NCurr(binW, :) ./ paddFourierFactor;
    [fxSCurr, fySCurr] = meshgrid((0:(Nt_SCurr(1) - 1)) / Nt_SCurr(1), (0:(Nt_SCurr(2) - 1)) / Nt_SCurr(2).');
    omSCurr = -2 * 1i * pi * ([fySCurr(:), fxSCurr(:)]);
    phaseSCurr = exp(omSCurr * nshiftSCurr(:)).';
    PHASE{binW}  = reshape(phaseSCurr, Nt_SCurr);
end
