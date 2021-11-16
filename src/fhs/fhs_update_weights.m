function [weights1, weights0] = fhs_update_weights(x_facet, size_v1, ...
    I, offset, status, nlevel, wavelet, Ncoefs, dims_overlap_ref, ...
    offsetL, offsetR, reweight_alpha, crop_sparsity, crop_low_rank, ...
    apodization_window, sig, sig_bar)
% Update the weights of the per facet priors.
%
% Update the weights involved in the reweighting procedure applied to the
% faceted l21 and nuclear norms. This version includes a spatial weighting
% correction for the nuclear norm (tapering window).
%
% Parameters
% ----------
% x_facet : array, double
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
% crop_sparsity : array, int
%     Relative cropping necessary for the facet l21-norm [1, 2].
% crop_low_rank : array, int
%     Relative cropping necessary for the facet nuclear norm [1, 2].
% apodization_window : array, double
%     Spatial weights applied to the facet nuclear norm (same size as
%     ``x_facet`` after cropping by ``crop_low_rank``).
% sig : double
%     Noise level (wavelet space).
% sig_bar : double
%     Noise level (singular value space).
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
% Last revised: [16/11/2021]
%-------------------------------------------------------------------------%
%%

%%  nuclear norm
sol = apodization_window .* x_facet(crop_low_rank(1) + 1:end, crop_low_rank(2) + 1:end, :);
sol = reshape(sol, [numel(sol) / size(sol, 3), size(x_facet, 3)]);
[~, z, ~] = svd(sol, 'econ'); clear sol;
z = abs(diag(z)); 
upsilon_bar = sig_bar * reweight_alpha;
weights0 = upsilon_bar ./ (upsilon_bar + z);  clear z;
%fprintf('\n nuclear weights done')


%% l21 norm
% ! offset for the zero-padding
spatial_size = dims_overlap_ref + offsetL + offsetR;
z = zeros(size_v1);
for l = 1:size(x_facet, 3)
    x_curr = zeros(spatial_size); % AD: mem reaons
    x_curr(offsetL(1) + 1:end - offsetR(1), offsetL(2) + 1:end - offsetR(2)) = x_facet(crop_sparsity(1) + 1:end, crop_sparsity(2) + 1:end, l);
    z(:, l) = sdwt2_sara_faceting(x_curr, I, offset, status, nlevel, wavelet, Ncoefs);
end, clear x_curr;
z = sqrt(sum(z.^2, 2)); 
upsilon = sig * reweight_alpha;
weights1 = upsilon ./ (upsilon + z); clear z;
%fprintf('\n l21  weights done')

end
