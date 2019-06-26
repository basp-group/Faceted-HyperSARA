function [v0, g0] = update_nuclear_spmd(v0, xhat, weights0, beta0)
%update_nuclear_spmd: update the dual variable realted to the nuclear norm
% prior.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > v0        dual variable associated with the nuclear nor prior [min(M*N, L), 1]
% > xhat      auxiliary variable related to the wideband image [M, N, L]
% > weights0  reweigthing weigths [min(M*N, L), 1]
% > beta0     thresholding parameter (1 / gamma0) [1]
%
% Output:
%
% < v0        dual variable associated with the nuclear nor prior [min(M*N, L), 1]
% < g0        auxiliary variable for the update of the primal variable [M, N, L]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [../../2019]
%-------------------------------------------------------------------------%
%%

[M, N, L] = size(xhat);
xhatm = reshape(xhat,numel(xhat)/L, L);
[U0, S0, V0] = svd(v0 + xhatm,'econ');
v0 = v0 + xhatm - (U0*diag(max(diag(S0) - beta0 * weights0, 0))*V0');
g0 = reshape(v0, M, N, L);

end