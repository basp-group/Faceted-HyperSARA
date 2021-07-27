function [mu_q, mu_bar_q] = compute_facet_log_prior_regularizer(x_overlap, Iq, ...
    offset, status_q, nlevel, wavelet, Ncoefs_q, dims_overlap_ref_q, ...
    offsetLq, offsetRq, crop_l21, crop_nuclear, w, size_v1, sig, sig_bar_q, regtype)
% Compute the value of the faceted prior (:math:`\ell_{2,1}` + nuclear 
% norm).
%
% Compute the value of the faceted prior (:math:`\ell_{2,1}` and  nuclear 
% norm) for a single facet. This version includes a spatial correction of 
% the faceted nuclear norm (tapering window).
%
% Parameters
% ----------
% x_overlap : array, double
%     Image facet, including overlap [M, N, L]
% Iq : array, int
%     Starting index of the non-overlapping facet [1, 2].
% offset : array, int
%     Offset to be used from one dictionary to another (different overlap
%     needed for each dictionary -> cropping) {nDictionaries}.
% status_q : cell, string
%     Status of the current facet (last or first facet along vert. / hrz.
%     direction).
% nlevel : int
%     Depth of the wavelet decompositions.
% wavelet : cell, string
%     Name of the wavelet dictionaries.
% Ncoefs_q : array, int
%     Size of the wavelet decompositions at each scale.
% dims_overlap_ref_q : array, int
%     Dimension of the facet [1, 2].
% offsetLq : array, int
%     Amount of zero-padding from the "left" [1, 2].
% offsetRq : array, int
%     Amount of zero-padding from the "right" [1, 2].
% crop_l21 : array, int
%     Relative cropping necessary for the facet :math:`\ell_[2,1]`-norm 
%     [1, 2].
% crop_nuclear : array, int
%     Relative cropping necessary for the facet nuclear norm [1, 2].
% w : array, double
%     Spatial weights applied to the facet nuclear norm (same size as 
%     ``x_overlap`` after cropping by ``crop_nuclear``).
% size_v1 : array, int
%     Size of the dual variable associated with the facet 
%     :math:`\ell_{2, 1}` norm [1, 2].
% sig : double
%     Regularization parameter for the :math:`\ell_{2, 1}` norm.
% sig_bar_q : double
%     Regularization parameter for the nuclear norm.
% regtype : string
%     [description]
%
% Returns
% -------
% mu_q : double
%     Regularization parameter for the :math:`\ell-{2, 1}` norm.
% mu_bar_q : double
%     Regularization parameter for the nuclear norm.
%
% .. note:: 
%
%   the l21 and nuclear norms do not act on the same facet, hence the
%   cropping step.
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [09/05/2021]
%-------------------------------------------------------------------------%
%% compute facet nuclear norm
c = size(x_overlap, 3);
xhatm = w.*x_overlap(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :);
xhatm = reshape(xhatm,numel(xhatm)/c,c);
[~,S0,~] = svd(xhatm,'econ');
% mu_bar_q = alph_bar / (sig_bar_q * sum(log(abs(diag(S0))/sig_bar_q + 1)));

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

switch regtype
    case "log"
        mu_bar_q = sig_bar_q * sum(log(abs(diag(S0))/sig_bar_q + 1));
        mu_q = sig * sum(log(sqrt(sum(abs(z).^2,2))/sig + 1));
    case "inv"
        mu_bar_q = sum(abs(diag(S0)));
        mu_q = sum(sqrt(sum(abs(z).^2,2)));
    otherwise
        mu_bar_q = sum(abs(diag(S0)));
        mu_q = sum(sqrt(sum(abs(z).^2,2)));
end

end