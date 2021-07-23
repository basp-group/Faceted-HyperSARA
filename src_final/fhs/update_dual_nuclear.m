function [v0, g0] = update_dual_nuclear(v0, xhat, w, weights0, beta0)
% Update the dual variable related to the weighted nuclear norm prior.
%
% Update step for the facet dual variable related to the nuclear norm 
% prior over a given facet. This version includes a spatial correction for 
% the faceted nuclear norm (tapering window).
%
% Args:
%     v0 (array): dual variable associated with the nuclear norm prior 
%                      [min(M*N, L), 1].
%     xhat (array): auxiliary variable related to the wideband image 
%                        [M, N, L].
%     w (array): spatial weights applied to the facet nuclear norm 
%                     [M, N, L].
%     weights0 (array): weights for the reweighting [min(M*N, L), 1].
%     beta0 (double): thresholding parameter (gamma0/sigma0). 
%
% Returns:
%     v0 (array): dual variable associated with the nuclear norm prior 
%                      [min(M*N, L), 1].
%     g0 (array): auxiliary variable for the update of the primal 
%                      variable [M, N, L].

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [08/08/2019]
%-------------------------------------------------------------------------%
%%

[M, N, c] = size(xhat);
xhatm = reshape(xhat.*w, numel(xhat)/c, c);
[U0, S0, V0] = svd(v0 + xhatm, 'econ');
v0 = v0 + xhatm - (U0*diag(max(diag(S0) - beta0 * weights0, 0))*V0');
g0 = w.*reshape(v0, M, N, c);

end
