function [weights1, weights0] = update_weights_overlap(x_overlap, size_v1, ...
    I, offset, status, nlevel, wavelet, Ncoefs, dims_overlap_ref, ...
    offsetL, offsetR, reweight_alpha, crop_l21, crop_nuclear, w)
% Update the weights of the per facet priors.
%
% Update the weights involved in the reweighting procedure applied to the 
% faceted l21 and nuclear norms. This version includes a spatial weighting 
% correction for the nuclear norm (tapering window).
%
% Parameters
% ----------
% x_overlap : array, double
%     Overlapping image facet [M, N, L].
% size_v1 : int
%     Size of the dual variable associated with the facet 
%     :math:`\ell_{2,1}` norm [1, 2].
% I : array,int
%     Starting index of the non-overlapping facet [1, 2].
% offset : array, int
%     Offset to be used from one dictionary to another (different overlap 
%     needed for each dictionary -> cropping) {nDictionaries}.
% status : array, int
%     Status of the current facet (last or first facet along vert. / hrz.
%     direction) [1, 2].
% nlevel : int
%     Depth of the wavelet decompositions.
% wavelet : cell, string
%     Name of the wavelet dictionaries.
% Ncoefs : array, int
%     Size of the wavelet decompositions at each scale.
% dims_overlap_ref : array, int
%     Dimension of the facet [1, 2].
% offsetL : array, int
%     Amount of zero-pading from the "left" [1, 2].
% offsetR : array, int
%     Amount of zero-padding from the "right" [1, 2].
% reweight_alpha : double
%     Reweighting parameter.
% crop_l21 : array, int
%     Relative cropping necessary for the facet l21-norm [1, 2].
% crop_nuclear : array, int
%     Relative cropping necessary for the facet nuclear norm [1, 2].
% w : array, double
%     Spatial weights applied to the facet nuclear norm (same size as 
%     ``x_overlap`` after cropping by ``crop_nuclear``).
%
% Returns
% -------
% weights1 : array, double
%      Weights for the reweighting of the facet :math:`\ell_{2,1}`-norm 
%      [size_v1(1), size_v1(2)].
% weights0 : array, double
%     Weights for the reweighting of the nuclear norm [min(M*N, L), 1]. 
%                    

%-------------------------------------------------------------------------%
%%
% Motivation: for loops are slow in spmd if not encapsulated in a function.
% Note: the l21 and nuclear norms do not act on the same facet, hence the
% cropping step.
%%
% Code: P.-A. Thouvenin.
% Last revised: [18/09/2020]
%-------------------------------------------------------------------------%
%%

% nuclear norm
sol = w.*x_overlap(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :);
sol = reshape(sol, [numel(sol)/size(sol, 3), size(x_overlap, 3)]);
[~,S00,~] = svd(sol,'econ');
d0 = abs(diag(S00));
% upsilon_bar = sig_bar*reweight_alpha;
% weights0 = upsilon_bar ./ (upsilon_bar + d0);
weights0 = reweight_alpha ./ (reweight_alpha + d0);

% l21 norm
zerosNum = dims_overlap_ref + offsetL + offsetR; % offset for the zero-padding
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetL(1)+1:end-offsetR(1), offsetL(2)+1:end-offsetR(2), :) = x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :);
z = zeros(size_v1);
for l = 1 : size(x_, 3)
    z(:, l) = sdwt2_sara_faceting(x_(:, :, l), I, offset, status, nlevel, wavelet, Ncoefs);
end
d1 = sqrt(sum(z.^2,2));
% upsilon = sig*reweight_alpha;
% weights1 = upsilon ./ (upsilon + d1);
weights1 = reweight_alpha ./ (reweight_alpha + d1);

end
