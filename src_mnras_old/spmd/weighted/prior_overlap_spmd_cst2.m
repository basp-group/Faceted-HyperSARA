function [l21_norm, nuclear_norm] = prior_overlap_spmd_cst2(x_overlap, Iq, ...
    dims_q, offset, status_q, nlevel, wavelet, Ncoefs_q, dims_overlap_ref_q, ...
    offsetLq, offsetRq, crop)
%prior_overlap_spmd_cst2: compute the value of the faceted prior (l21 + 
% nuclear norm) for a single facet, assuming the faceted nuclear norm acts
% on a larger overlap than faceted wavelet transforms.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x_overlap               overlapping image facet [M, N, L]
% > Iq                      starting index of the non-overlapping facet 
%                           [1, 2]
% > dims_q                  size of the non-overlapping facet [1, 2] 
% > offset                  offset to be used from one dictionary to
%                           another (different overlap needed for each 
%                           dictionary -> cropping) {nDictionaries}
% > status_q                status of the current facet (vert. and hrz) ...
% > nlevel                  depth of the wavelet decompositions
% > wavelet                 name of the wavelet dictionaries
% > Ncoefs_q                size of the wavelet decompositions at each
%                           scale
% > dims_overlap_ref_q      size of the facet [1, 2] 
% > offsetLq                amount of zero-pading from the "left" [1, 2]
% > offsetRq                amout of zero-padding from the "right" [1, 2]
% > crop                    relative cropping necessary for the nuclear 
%                           norm (smaller overlap than the l21) [1, 2]
% 
% Note: the l21 and nuclear norms do not act on the same facet, hence the
% cropping step.
%
% Output:
%
% < l21_norm                facet l21-norm [1]
% < nuclear_norm            facet nuclear norm [1]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%% compute facet nuclear norm
c = size(x_overlap, 3);
xhatm = reshape(x_overlap,numel(x_overlap)/c,c);
[~,S0,~] = svd(xhatm,'econ');
nuclear_norm = norm(diag(S0),1);

%% compute facet l21-norm
% zero-padding
zerosNum = dims_overlap_ref_q + offsetLq + offsetRq;
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap(crop(1)+1:end, crop(2)+1:end, :);

% faceted SARA
l21_norm = 0;
for l = 1 : c
    wl = sdwt2_sara(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
    l2 = sqrt(sum(abs(wl).^2,2));
    l21_norm = l21_norm + norm(l2(:),1);
end

end
