function [v1, weights1, s] = hs_initialize_dual_sparsity_distributed(x, ...
    Psit, extension_mode, nlevel, reweighting_alpha, sig)
% Initalize the dual variables related to the
% l21-norm prior.
%
% Parameters
% ----------
% x : array (3d)
%     Wideband image [M, N, L].
% Psit : anonymous function
%     Full SARA operator @[1].
% extension_mode : string
%     Name of the boundary extension mode.
% nlevel : int
%     Depth of the wavelet decompositions.
% reweighting_alpha : double
%     Rewighting parameter.
% sig : double
%     Estimate of the noise level transferred to the SARA domain.
%
% Returns
% -------
% v1 : array (double, 2d)
%     Dual variable associated with the l21-norm prior [s, L].
% weights1 : array (double, 2d)
%     Weights associated for the reweigthing step [s, L].
% s : int
%     Number of wavelet decompostion for each channel [1].
%

% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
% -------------------------------------------------------------------------%
%%

[M, N, c] = size(x);
% number of cofficients resulting from the 8 first Daubcehies wavelet
% transforms
[~, s] = n_wavelet_coefficients(2 * (1:8)', [M, N], extension_mode, nlevel);
% add size from Dirac dictionary
s = s + M * N;

v1 = zeros(s, c);
for l = 1:c
    v1(:, l) = Psit(x(:, :, l));
end
% compute part of the row-wise l2 norm
d1_ = sum(v1.^2, 2);
% compute sum across all workers
d1 = gplus(d1_);
d1 = sqrt(d1);
upsilon = sig * reweighting_alpha;
weights1 = upsilon ./ (upsilon + d1);

v1 = zeros(s, c);

end
