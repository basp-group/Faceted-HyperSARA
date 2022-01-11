function [l21_norm, nuclear_norm] = ...
    fhs_compute_facet_prior(x_facet, Iq, offset, status_q, ...
    nlevel, wavelet, Ncoefs_q, dims_overlap_ref_q, offsetLq, offsetRq, ...
    crop_sparsity, crop_low_rank, spatial_weights, size_v1)
% Compute the value of the faceted prior (:math:`\ell_{2,1}` + nuclear
% norm).
%
% Compute the value of the faceted prior (:math:`\ell_{2,1}` and  nuclear
% norm) for a single facet. This version includes a spatial correction of
% the faceted nuclear norm (tapering window).
%
% Parameters
% ----------
% x_facet : double[:, :, :]
%     Overlapping image facet ``[M, N, L]``.
% Iq : int[:]
%     Starting index of the non-overlapping facet ``[1, 2]``.
% offset : int[:]
%     offset to be used from one dictionary to another (different overlap
%     needed for each dictionary -> cropping) ``{nDictionaries}``.
% status_q : int[:]
%     Status of the current facet (last or first facet along vert. / hrz.
%     direction).
% nlevel : int
%     Depth of the wavelet decompositions.
% wavelet : cell (strings)
%     Name of the wavelet dictionaries.
% Ncoefs_q : int[:, :]
%     Size of the wavelet decompositions at each scale.
% dims_overlap_ref_q : int[:]
%     Dimensions of the facet ``[1, 2]``.
% offsetLq : int[:]
%     Amount of zero-padding from the "left" ``[1, 2]``.
% offsetRq : int[:]
%     Amount of zero-padding from the "right" ``[1, 2]``.
% crop_sparsity : int[:]
%     Relative cropping necessary for the facet :math:`\ell_{2,1}` norm
%     [1, 2].
% crop_low_rank : int[:]
%     Relative cropping necessary for the facet nuclear norm ``[1, 2]``.
% spatial_weights : double[:, :]
%     Spatial weights (apodization window) applied to the facet nuclear
%     norm.
% size_v1 : int[:]
%     Size of the dual variable associated with the facet
%     :math:`\ell_{2,1}` norm ``[1, 2]``.
%
% Returns
% -------
% l21_norm  : double
%     Facet :math:`\ell_{2,1}` norm.
% nuclear_norm  : double
%     Facet nuclear norm.
%
% ..note::
%
%   The l21 and nuclear norms do not act on the same facet, hence the
%   cropping step.
%

% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
% -------------------------------------------------------------------------%
%% compute facet nuclear norm
n_channels = size(x_facet, 3);
xhatm = spatial_weights .* x_facet(crop_low_rank(1) + 1:end, crop_low_rank(2) + 1:end, :);
xhatm = reshape(xhatm, numel(xhatm) / n_channels, n_channels);
[~, S0, ~] = svd(xhatm, 'econ');      clear xhatm;
nuclear_norm = norm(diag(S0), 1);     clear S0;

%% compute facet l21-norm
% zero-padding
zerosNum = dims_overlap_ref_q + offsetLq + offsetRq;

% faceted SARA
z = zeros(size_v1);
for l = 1:n_channels
    x_curr = zeros(zerosNum);
    x_curr(offsetLq(1) + 1:end - offsetRq(1), offsetLq(2) + 1:end - offsetRq(2)) =  x_facet(crop_sparsity(1) + 1:end, crop_sparsity(2) + 1:end, l);
    z(:, l) = sdwt2_sara_faceting(x_curr, Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
end; clear x_curr;
l21_norm = sum(sqrt(sum(abs(z).^2, 2)), 1); clear z;

end
