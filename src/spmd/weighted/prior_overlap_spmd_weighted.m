function [l21_norm, nuclear_norm] = prior_overlap_spmd_weighted(x_overlap, ...
    Iq, offset, status_q, nlevel, wavelet, Ncoefs_q, dims_overlap_ref_q, ...
    offsetLq, offsetRq, w, size_v1)
%prior_overlap_spmd_weighted: compute the value of the faceted prior 
% (l21 + nuclear norm) for a single facet. This version includes a
% spatial correction for the faceted nuclear norm (tapering window).
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x_overlap               overlapping image facet [M, N, L]
% > Iq                      starting index of the non-overlapping facet 
%                           [1, 2]
% > offset                  offset to be used from one dictionary to
%                           another (different overlap needed for each 
%                           dictionary -> cropping) {nDictionaries}
% > status_q                status of the current facet (last or first 
%                           facet along vert. or hrz. direction)
% > nlevel                  depth of the wavelet decompositions
% > wavelet                 name of the wavelet dictionaries
% > Ncoefs_q                size of the wavelet decompositions at each
%                           scale
% > dims_overlap_ref_q      dimension of the facet [1, 2] 
% > offsetLq                amount of zero-pading from the "left" [1, 2]
% > offsetRq                amout of zero-padding from the "right" [1, 2]
%                           norm [size(x_overlap)] 
% > w                       spatial weights applied to the facet nuclear
%                           norm (same size as x_overlap)
% > size_v1                 size of the dual variable associated with the 
%                           facet l21-norm [1, 2]
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
xhatm = reshape(x_overlap.*w,prod(dims_overlap_ref_q),c);
[~,S0,~] = svd(xhatm,'econ');
nuclear_norm = norm(diag(S0),1);

%% compute facet l21-norm
% zero-padding
zerosNum = dims_overlap_ref_q + offsetLq + offsetRq;
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap;

% faceted SARA
z = zeros(size_v1);
for l = 1 : c
    z(:,l) = sdwt2_sara(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
end
l21_norm = sum(sqrt(sum(abs(z).^2,2)), 1);

end
