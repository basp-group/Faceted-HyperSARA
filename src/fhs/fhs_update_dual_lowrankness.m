function [v0, g0] = fhs_update_dual_lowrankness(v0, xhat, ...
    apodization_window, weights0, beta0)
% Update the dual variable related to the weighted nuclear norm prior.
%
% Update step for the facet dual variable related to the nuclear norm
% prior over a given facet. This version includes a spatial correction for
% the faceted nuclear norm (tapering window).
%
% Parameters
% ----------
% v0 : double[:]
%     Dual variable associated with the nuclear norm prior
%     ``[min(M*N, L), 1]``.
% xhat : double[:, :, :]
%     Auxiliary variable related to the wideband image ``[M, N, L]``.
% apodization_window : double[:, :]
%     Spatial weights applied to the facet nuclear norm ``[M, N, L]``.
% weights0 : double[:]
%     Weights for the reweighting ``[min(M*N, L), 1]``.
% beta0 : double
%     Thresholding parameter (gamma0/sigma0).
%
% Returns
% -------
% v0 : double[:, :]
%     Dual variable associated with the nuclear norm prior
%     ``[min(M*N, L), 1]``.
% g0 : double[:, :, :]
%     Auxiliary variable for the update of the primal variable
%     ``[M, N, L]``.
%

% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [08/08/2019]
% -------------------------------------------------------------------------%
%%

[M, N, c] = size(xhat);
xhat_plus_v0 = v0 + reshape(xhat .* apodization_window, numel(xhat) / c, c);
[U0, S0, V0] = svd(xhat_plus_v0, 'econ');
S0 = diag(S0);
v0 = xhat_plus_v0  - (U0 * diag(max(S0 - beta0 * weights0, 0)) * V0');  clear S0 U0 V0 xhat_plus_v0;
g0 = apodization_window .* reshape(v0, M, N, c);

end
