function [l21_norm, nuclear_norm] = ...
    fhs_compute_facet_prior(x_overlap, Iq, offset, status_q, ...
    nlevel, wavelet, Ncoefs_q, dims_overlap_ref_q, offsetLq, offsetRq, ...
    crop_l21, crop_nuclear, spatial_weights, size_v1)    
% Compute the value of the faceted prior (:math:`\ell_{2,1}` + nuclear 
% norm).
%
% Compute the value of the faceted prior (:math:`\ell_{2,1}` and  nuclear 
% norm) for a single facet. This version includes a spatial correction of 
% the faceted nuclear norm (tapering window).
%
% Parameters
% ----------
% x_overlap : array (3d)
%     Overlapping image facet [M, N, L].
% Iq : array (2d)
%     Starting index of the non-overlapping facet [1, 2].
% offset : array (1d)
%     offset to be used from one dictionary to another (different overlap
%     needed for each dictionary -> cropping) {nDictionaries}.
% status_q : array (1d)
%     Status of the current facet (last or first facet along vert. / hrz.
%     direction).
% nlevel : int
%     Depth of the wavelet decompositions.
% wavelet : cell (strings)
%     Name of the wavelet dictionaries.
% Ncoefs_q : array (2d)
%     Size of the wavelet decompositions at each scale.
% dims_overlap_ref_q : array (1d)
%     Dimensions of the facet [1, 2].
% offsetLq : array (1d)
%     Amount of zero-padding from the "left" [1, 2].
% offsetRq : array (1d)
%     Amount of zero-padding from the "right" [1, 2].
% crop_l21 : array (1d)
%     Relative cropping necessary for the facet :math:`\ell_{2,1}` norm 
%     [1, 2].
% crop_nuclear : array (1d)
%     Relative cropping necessary for the facet nuclear norm [1, 2].
% spatial_weights : array (2d)
%     Spatial weights (apodization window) applied to the facet nuclear 
%     norm.
% size_v1 : array (1d)
%     Size of the dual variable associated with the facet 
%     :math:`\ell_{2,1}` norm [1, 2].
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

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%% compute facet nuclear norm
c = size(x_overlap, 3);
xhatm = spatial_weights.*x_overlap(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :);
xhatm = reshape(xhatm,numel(xhatm)/c,c);
[~,S0,~] = svd(xhatm,'econ');
nuclear_norm = norm(diag(S0),1);

%% compute facet l21-norm
% zero-padding
zerosNum = dims_overlap_ref_q + offsetLq + offsetRq;
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :);

% faceted SARA
z = zeros(size_v1);
for l = 1 : c
    z(:,l) = sdwt2_sara_faceting(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
end
l21_norm = sum(sqrt(sum(abs(z).^2,2)), 1);

end
