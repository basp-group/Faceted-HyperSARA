function [v1, g1] = update_dual_l21_serial(v1, Psit, Psi, xhat, weights1, beta1, sigma1)
%update_dual_nuclear_serial: update the dual variable related to the facet l21-norm
% prior.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > v1                      dual variable associated with the facet l21-
%                           norm [s, L]
% > Psit                    full SARA operator @[1]                     
% > Psi                     adjoint of the full SARA operator @[1] 
% > xhat                    auxiliary variable related to the wideband
%                           image [M, N, L]
% > weights1                weights associated with the reweigthing step
%                           pixels due to the overlap between facets [M, N]
% > beta1                   update step (mu / gamma1) [1]
% > sigma1                  convergence parameter [1]
%
% Output:
%
% < v1                      dual variable associated with the l21-norm 
%                           prior [s, L]
% < g1                      auxiliary variable for the update of the primal
%                           variable [M, N, L]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%

r1 = zeros(size(v1));
g1 = zeros(size(xhat));

for l = 1:size(xhat, 3)
    r1(:, l) = v1(:, l) + Psit(xhat(:,:,l));
end
l2 = sqrt(sum(abs(r1).^2,2));
l2_soft = max(l2 - beta1*weights1, 0)./ (l2+eps);
v1 = r1 - (l2_soft .* r1);

for l = 1:size(xhat, 3)
    g1(:,:,l) = sigma1*Psi(v1(:,l));
end

end
