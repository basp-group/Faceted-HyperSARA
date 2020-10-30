function [v0, g0] = update_nuclear_spmd(v0, xhat, weights0, beta0)
% Update the dual variable related to the nuclear norm prior.
%
% Update step for the facet dual variable related to the nuclear norm 
% prior over a given facet.
%
% Args:
%     v0 (array_like): dual variable associated with the nuclear norm prior 
%                      [min(M*N, L), 1].
%     xhat (array_like): auxiliary variable related to the wideband image 
%                        [M, N, L],
%     weights0 (array_like): weights for the reweighting [min(M*N, L), 1].
%     beta0 (double): thresholding parameter (gamma0/sigma0). 
%
%
% Returns:
%     v0 (array_like): dual variable associated with the nuclear norm prior 
%                      [min(M*N, L), 1].
%     g0 (array_like): auxiliary variable for the update of the primal 
%                      variable [M, N, L].

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%

[M, N, L] = size(xhat);
xhatm = reshape(xhat,numel(xhat)/L, L);
[U0, S0, V0] = svd(v0 + xhatm,'econ');
v0 = v0 + xhatm - (U0*diag(max(diag(S0) - beta0 * weights0, 0))*V0');
g0 = reshape(v0, M, N, L);

end
