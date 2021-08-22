function [xsol, xhat, rel_x, norm_x] = update_primal(xsol, g)
% Update the wideband facet (primal variable).
%
% Update step for the primal variable in the preconditioned primal-dual 
% algorithm.
%
% Parameters
% ----------
% xsol : array
%     Wideband facet [M, N, L].
% g : array
%     Auxiliary variable related to the wideband facet [M, N, L].
%
% Returns
% -------
% xsol : array
%     Udated wideband facet [M, N, L].
% xhat : array
%     Auxiliary variable related to the wideband facet [M, N, L].
% rel_x : double
%     Squared relative variation b etween two consecutive iterates.
% norm_x : double
%     Squared Euclidean norm of `xsol`.
%

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%

prev_xsol = xsol;
xsol = max(real(xsol - g), 0);
xhat = 2*xsol - prev_xsol;

rel_x = sum(power(prev_xsol(:) - xsol(:), 2));
norm_x = sum(power(xsol(:), 2));

end
