function [v0, g0] = hs_update_dual_lowrankness_serial(v0, xhat, ...
    weights0,  beta0, sigma0)
% Update the dual variables related to the nuclear
% norm prior.
%
% Parameters
% ----------
% v0 : double[:, :]
%     Dual variable associated with the nuclear norm prior
%     ``[min(M*N, L), 1]``.
% xhat : double[:, :, :]
%     Auxiliary variable related to the wideband image ``[M, N, L]``.
% weights0 : double[:, :]
%     Weights for the reweighting ``[min(M*N, L), 1]``.
% beta0 : double
%     Thresholding parameter (gamma0 / sigma0) [1]
% sigma0 : double
%     Convergence parameter [1].
%
% Returns
% -------
% v0 : double[:, :]
%     Dual variable associated with the nuclear norm prior
%     ``[min(M*N, L), 1]``.
% weights0 : double[:, :, :]
%     Auxiliary variable for the update of the primal variable
%     ``[M, N, L].``
%

% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
% -------------------------------------------------------------------------%
%%

[M, N, L] = size(xhat);
xhatm = reshape(xhat, M * N, L);
[U0, S0, V0] = svd(v0 + xhatm, 'econ');
v0 = v0 + xhatm - (U0 * diag(max(diag(S0) - beta0 * weights0, 0)) * V0');
clear V0 U0 S0 xhatm;
g0 = sigma0 * reshape(v0, M, N, L);

end
