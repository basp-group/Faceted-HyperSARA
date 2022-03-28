function [v1, g1] = hs_update_dual_sparsity_distributed(v1, Psit, Psi, ...
    xhat, weights1, beta1, sigma1)
% Update the dual variable related to the facet :math:`\ell_{2,1}` norm
% prior.
%
% Parameters
% ----------
% v1 : double[:, :]
%     Dual variable associated with the facet :math:`\ell_{2,1}` norm
%     ``[s, L]``.
% Psit : anonymous function
%     Full SARA operator ``@[1]``.
% Psi : anonymous function
%     Adjoint of the full SARA operator ``@[1]``.
% xhat : double[:, :, :]
%     Auxiliary variable related to the wideband image ``[M, N, L]``.
% weights1 : double[:, :]
%     Weights associated with the reweigthing step pixels due to the
%     overlap between facets ``[M, N]``.
% beta1 : double
%     Update step (mu / gamma1) ``[1]``.
% sigma1 : double
%     Convergence parameter ``[1]``.
%
% Returns
% -------
% v1 : double[:, :]
%     Dual variable associated with the :math:`\ell_{2,1}` norm prior
%     ``[s, L]``.
% g1 : double[:, :, :]
%     Auxiliary variable for the update of the primal variable
%     ``[M, N, L]``.
%

% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
% -------------------------------------------------------------------------%
%%

r1 = zeros(size(v1));
g1 = zeros(size(xhat));

for l = 1:size(xhat, 3)
    r1(:, l) = v1(:, l) + Psit(xhat(:, :, l));
end
l2_ = sum(abs(r1).^2, 2);
l2 = gplus(l2_);
l2 = sqrt(l2);
l2_soft = max(l2 - beta1 * weights1, 0) ./ (l2 + eps);
v1 = r1 - (l2_soft .* r1);

for l = 1:size(xhat, 3)
    g1(:, :, l) = sigma1 * Psi(v1(:, l));
end

end
