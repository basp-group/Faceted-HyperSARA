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
% v1 : double[:, :]
%     Dual variable associated with the :math:`\ell_{2,1}` prior ``[s, L]``.
% x_facet : double[:, :, :]
%     Overlapping image facet ``[M, N, L]``.
% apodization_window : double[:, :]
%     Weights to mitigate tessellation aretefacts due to the overlap
%     between facets ``[M, N]``.
% beta1 : double
%     Ratio between regularization and convergence parameter
%     (gamma1 / sigma1).
% Iq : int[:]
%     Starting index of the non-overlapping base facet ``[1, 2]``.
% dims_q : int[:]
%     Dimensions of the non-overlapping base facet ``[1, 2]``.
% I_overlap_q : int[:]
%     Starting index of the facet ``[1, 2]``.
% dims_overlap_q : int[:]
%     Dimensions of the facet ``[1, 2]``.
% offset : cell of int[:]
%     Offset to be used from one dictionary to another (different overlap
%     needed for each dictionary -> cropping) ``{nDictionaries}``.
% status_q : int[:]
%     Status of the current facet (last or first facet along vert. or hrz.
%     direction) [ndict, 2].
% nlevel : int
%     Depth of the wavelet decompositions.
% wavelet : cell (string)
%     Name of the wavelet dictionaries {ndict}.
% Ncoefs_q : int[:]
%     Size of the wavelet decompositions at each scale.
% temLIdxs_q : int[:]
%     Amount of cropping from the "left" ``[1, 2]``.
% temRIdxs_q : int[:]
%     Amount of cropping from the "right" ``[1, 2]``.
% offsetLq : int[:]
%     Amount of zero-pading from the "left" ``[1, 2]``.
% offsetRq : int[:]
%     Amount of zero-pading from the "right" ``[1, 2]``.
% dims_overlap_ref_q : int[:]
%     Dimension of the facet ``[1, 2]``.
%
% Returns
% -------
% v1 : double[:, :]
%     Dual variable associated with the :math:`\ell_{2,1}` prior ``[s, L]``.
% g : double[:, :, :]
%     Auxiliary variable for the update of the primal variable ``[M, N, L]``.
%

% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revision: [30/11/2020]
% -------------------------------------------------------------------------%
%%

spatial_size = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding
spectral_size = size(x_facet, 3);
x_facet_size = size(x_facet);

x_curr = zeros(spatial_size);
x_curr(offsetLq(1) + 1:end - offsetRq(1), offsetLq(2) + 1:end - offsetRq(2)) = x_facet(:, :, 1);
x_curr = sdwt2_sara_faceting(x_curr, Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);

w = zeros(numel(x_curr), spectral_size);
w(:, 1) = x_curr; clear x_curr;
for l = 2:spectral_size
    x_curr = zeros(spatial_size);
    x_curr(offsetLq(1) + 1:end - offsetRq(1), offsetLq(2) + 1:end - offsetRq(2)) = x_facet(:, :, l);
    w(:, l) = sdwt2_sara_faceting(x_curr, Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
end; clear  x_curr;

w = v1 +  w;

l2 = sqrt(sum(abs(w).^2, 2));
l2_soft = max(l2 - beta1 * apodization_window, 0) ./ (l2 + eps);

v1 = w - l2_soft .* w;  clear w;

g = zeros(x_facet_size);
for l = 1:spectral_size
    g(:, :, l) = isdwt2_sara_faceting(v1(:, l), Iq, dims_q, I_overlap_q, dims_overlap_q, Ncoefs_q, nlevel, wavelet, temLIdxs_q, temRIdxs_q);
end

end
