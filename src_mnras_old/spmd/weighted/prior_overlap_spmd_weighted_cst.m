function [l21_norm, nuclear_norm] = prior_overlap_spmd_weighted_cst(x_overlap, Iq, ...
    dims_q, offset, status_q, nlevel, wavelet, Ncoefs_q, dims_overlap_ref_q, ...
    offsetLq, offsetRq, w, crop)
%prior_overlap_spmd_weighted_cst: compute the value of the faceted prior 
% (l21 + nuclear norm) for a single facet. This version includes a
% spatial correction for the faceted nuclear norm (tapering window).
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
% > status_q                status of the current facet (last or first 
%                           facet along vert. or hrz. direction)
% > nlevel                  depth of the wavelet decompositions
% > wavelet                 name of the wavelet dictionaries
% > Ncoefs_q                size of the wavelet decompositions at each
%                           scale
% > dims_overlap_ref_q      size of the facet [1, 2] 
% > offsetLq                amount of zero-pading from the "left" [1, 2]
% > offsetRq                amout of zero-padding from the "right" [1, 2]
%                           norm [size(x_overlap)] 
% > crop                    relative cropping necessary for the nuclear 
%                           norm (overlap for l21 is larger) [1, 2]
% > w                       spatial weights applied to the facet nuclear
%                           norm (same size as x_overlap after cropping by 
%                           crop)
% > size_v1                 size of the dual variable associated with the 
%                           facet l21-norm [1, 2]
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
xhatm = x_overlap(crop(1)+1:end, crop(2)+1:end, :).*w;
xhatm = reshape(xhatm,numel(xhatm)/c,c);
[~,S0,~] = svd(xhatm,'econ');
nuclear_norm = norm(diag(S0),1);

%% compute facet l21-norm
% zero-padding
zerosNum = dims_overlap_ref_q + offsetLq + offsetRq;
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap;

% faceted SARA
l21_norm = 0;
for l = 1 : c
    wl = sdwt2_sara(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
    l2 = sqrt(sum(abs(wl).^2,2));
    l21_norm = l21_norm + norm(l2(:),1);
end

end
