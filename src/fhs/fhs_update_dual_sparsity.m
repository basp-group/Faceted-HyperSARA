function [v1, g] = fhs_update_dual_sparsity(v1, x_facet, ...
    apodization_window, beta1, Iq, dims_q, I_overlap_q, ...
    dims_overlap_q, offset, status_q, nlevel, wavelet, Ncoefs_q, ...
    temLIdxs_q, temRIdxs_q, offsetLq, offsetRq, dims_overlap_ref_q)
% Update the dual variable associated with the facet l21-norm prior.
%
% Update a facet dual variable associated with the l21-norm prior.
%
% Parameters
% ----------
% v1 : array (2d)
%     Dual variable associated with the :math:`\ell_{2,1}` prior [s, L].
% x_facet : array (3d)
%     Overlapping image facet [M, N, L].
% apodization_window : array (2d)
%     Weights to mitigate tessellation aretefacts due to the overlap
%     between facets [M, N].
% beta1 : double
%     Ratio between regularization and convergence parameter
%     (gamma1 / sigma1).
% Iq : array (1d)
%     Starting index of the non-overlapping base facet [1, 2].
% dims_q : array (1d)
%     Dimensions of the non-overlapping base facet [1, 2].
% I_overlap_q : array (1d)
%     Starting index of the facet [1, 2].
% dims_overlap_q : array (1d)
%     Dimensions of the facet [1, 2].
% offset : cell
%     Offset to be used from one dictionary to another (different overlap
%     needed for each dictionary -> cropping) {nDictionaries}.
% status_q : array (1d)
%     Status of the current facet (last or first facet along vert. or hrz.
%     direction) [ndict, 2].
% nlevel : int
%     Depth of the wavelet decompositions.
% wavelet : cell (string)
%     Name of the wavelet dictionaries {ndict}.
% Ncoefs_q : array (1d)
%     Size of the wavelet decompositions at each scale.
% temLIdxs_q : array (1d)
%     Amount of cropping from the "left" [1, 2].
% temRIdxs_q : array (1d)
%     Amount of cropping from the "right" [1, 2].
% offsetLq : array (1d)
%     Amount of zero-pading from the "left" [1, 2].
% offsetRq : array (1d)
%     Amount of zero-pading from the "right" [1, 2].
% dims_overlap_ref_q : array (1d)
%     Dimension of the facet [1, 2].
%
% Returns
% -------
% v1 : array (2d)
%     Dual variable associated with the :math:`\ell_{2,1}` prior [s, L].
% g : array (3d)
%     Auxiliary variable for the update of the primal variable [M, N, L].
%

% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revision: [30/11/2020]
% -------------------------------------------------------------------------%
%%

spatial_size = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding
x_ = zeros([spatial_size, size(x_facet, 3)]);
x_(offsetLq(1) + 1:end - offsetRq(1), offsetLq(2) + 1:end - offsetRq(2), :) = x_facet;
g = zeros(size(x_facet));

tmp = sdwt2_sara_faceting(x_(:, :, 1), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
w = zeros(numel(tmp), size(x_, 3));
w(:, 1) = tmp;
for l = 2:size(x_, 3)
    w(:, l) = sdwt2_sara_faceting(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
end
w = v1 +  w;
l2 = sqrt(sum(abs(w).^2, 2));
l2_soft = max(l2 - beta1 * apodization_window, 0) ./ (l2 + eps);
v1 = w - l2_soft .* w;

for l = 1:size(x_, 3)
    g(:, :, l) = isdwt2_sara_faceting(v1(:, l), Iq, dims_q, I_overlap_q, dims_overlap_q, Ncoefs_q, nlevel, wavelet, temLIdxs_q, temRIdxs_q);
end

end
