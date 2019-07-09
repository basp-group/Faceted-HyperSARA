function [v0, g0] = update_nuclear_spmd_weighted(v0, xhat, w, weights0, beta0)
%update_nuclear_spmd: update the dual variable realted to the nuclear norm
% prior. This version includes the spatial correction of the faceted 
% nuclear norm.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > v0        dual variable associated with the nuclear norm prior 
%             [min(M*N, L), 1]
% > xhat      auxiliary variable related to the wideband image [M, N, L]
% > w         spatial weights applied to the facet nuclear norm (same size 
%             as xhat)
% > weights0  weights for the reweighting [min(M*N, L), 1]
% > beta0     thresholding parameter (1 / gamma0) [1]
%
% Output:
%
% < v0        dual variable associated with the nuclear norm prior 
%             [min(M*N, L), 1]
% < g0        auxiliary variable for the update of the primal variable 
%             [M, N, L]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [../../2019]
%-------------------------------------------------------------------------%
%%

[M, N, c] = size(xhat);
xhatm = reshape(xhat.*w, numel(xhat)/c, c);
[U0, S0, V0] = svd(v0 + xhatm, 'econ');
v0 = v0 + xhatm - (U0*diag(max(diag(S0) - beta0 * weights0, 0))*V0');
g0 = w.*reshape(v0, M, N, c);

end
