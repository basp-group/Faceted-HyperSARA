function [v0, v1, weights0, weights1] = ...
    fhs_initialize_dual_and_weights(x_facet, I, offset, status, nlevel, ...
    wavelet, Ncoefs, dims_o, n_channels, dims_overlap_ref, offsetL, offsetR, ...
    reweight_alpha, crop_sparsity, crop_low_rank, apodization_window, sig, sig_bar)
% Initialize dual variables and weights associated with the faceted
% low-rankness and average joint-sparsity prior.
%
% Parameters
% ----------
% x_facet : array, double
%     Overlapping image facet [M, N, L].
% I : array, int
%     Starting index of the non-overlapping facet [1, 2].
% offset : array, double
%     Offset to be used from one dictionary to another (different overlap
%     needed for each dictionary -> cropping) {nDictionaries}.
% status : array, int
%     Status of the current facet (last or first facet along vert. / hrz.
%     direction) [1, 2].
% nlevel : int
%     Depth of the wavelet decompositions.
% wavelet : cell, string
%     Name of the wavelet transforms involved in the faceted average
%     joint-sparsity prior.
% Ncoefs : array, double
%     Size of the wavelet decompositions at each scale.
% dims_o : array, double
%     Dimension of a facet (with overlap) [1, 2].
% n_channels : int
%     Number of spectral channels.
% dims_overlap_ref : array, int
%     Dimension of the facet [1, 2].
% offsetL : array, int
%     Amount of zero-pading from the "left" [1, 2].
% offsetR : array, int
%     Amount of zero-padding from the "right" [1, 2].
% reweight_alpha : double
%     Reweighting parameter.
% crop_sparsity : array, int
%     Amount of pixels to be cropped from the facet along each dimension to
%     retrieve the pixels over which the current facet's sparstiy prior is
%     acting [1, 2].
% crop_low_rank : array, int
%     Amount of pixels to be cropped from the facet along each dimension to
%     retrieve the pixels over which the current facet's sparstiy prior is
%     acting [1, 2].
% apodization_window : array, double
%     Apodization window used in the faceted low-rankness prior.
% sig : double
%     Noise level for the weights (joint-sparsity prior).
% sig_bar : double
%     Noise level for the weights (low-rankness prior).
%
% Returns
% -------
% v0 : array, double
%     Dual variable associated with the low-rankness prior.
% v1 : array, double
%     Dual variable associated with the low-rankness prior.
% weights0 : array, double
%     Weigths associated with the low-rankness prior.
% weights1 : array, double
%     Weigths associated with the average joint-sparsity prior.
%

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [../../....]
%-------------------------------------------------------------------------%
%%

% ! dual variables initialized to 0, and weights initialized from the
% ! current value of the primal variable

% ! checking if x_facet = 0 is not strictly speaking (weights are automatically
% ! 1 in this case), but it can save some operations

% compute size of the dual variables and the weigths
p = prod(Ncoefs, 2);
% if dirac_present
%     sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) + prod(dims);
% else
%     sz = 3*sum(p) - 2*sum(p(nlevel+1:nlevel+1:end));
% end
% number of wavelet coeffs when the Dirac basis is used
sz = 3 * sum(p(1:end)) - 2 * sum(p(nlevel + 1:nlevel + 1:end)) - 2 * p(end);
flag_zero = (norm(x_facet(:)) == 0);

% nuclear norm
v0 = zeros(prod(dims_o), n_channels);
if flag_zero
    weights0 = ones(min(prod(dims_o), n_channels), 1);
else
    sol = apodization_window .* x_facet(crop_low_rank(1) + 1:end, crop_low_rank(2) + 1:end, :);
    sol = reshape(sol, [numel(sol) / size(sol, 3), size(x_facet, 3)]);
    [~, S00, ~] = svd(sol, 'econ');
    d0 = abs(diag(S00));
    upsilon_bar = sig_bar * reweight_alpha;
    weights0 = upsilon_bar ./ (upsilon_bar + d0);
end

% l21 norm
v1 = zeros(sz, n_channels);
if flag_zero
    weights1 = ones(sz, 1);
else
    zerosNum = dims_overlap_ref + offsetL + offsetR; % offset for the zero-padding
    x_ = zeros([zerosNum, size(x_facet, 3)]);
    x_(offsetL(1) + 1:end - offsetR(1), offsetL(2) + 1:end - offsetR(2), :) = x_facet(crop_sparsity(1) + 1:end, crop_sparsity(2) + 1:end, :);
    z = zeros(sz, n_channels);
    for l = 1:size(x_, 3)
        z(:, l) = sdwt2_sara_faceting(x_(:, :, l), I, offset, status, nlevel, wavelet, Ncoefs);
    end
    d1 = sqrt(sum(z.^2, 2));
    upsilon = sig * reweight_alpha;
    weights1 = upsilon ./ (upsilon + d1);
end

end
