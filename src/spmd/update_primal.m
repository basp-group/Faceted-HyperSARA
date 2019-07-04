function [xsol, xhat, rel_x, norm_x] = update_primal(xsol, g)
%update_primal: update wideband image (primal variable).
%-------------------------------------------------------------------------%
%%
% Input:
%
% > xsol      wideband facet [M, N, L]
% > g         auxiliary variable related to the wideband facet [M, N, L]
%
% Output:
%
% < xsol      wideband facet [M, N, L]
% < xhat      auxiliary variable related to the wideband facet [M, N, L]
% < rel_x     square of the relative variation [1]
% < norm_x    square Euclidean norm of xsol [1] 
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [../../2019]
%-------------------------------------------------------------------------%
%%

prev_xsol = xsol;
xsol = max(real(xsol - g), 0);
xhat = 2*xsol - prev_xsol;

rel_x = sum(power(prev_xsol(:) - xsol(:), 2));
norm_x = sum(power(xsol(:), 2));

end
