function [xsol, xhat, rel_x, norm_x] = update_primal(xsol, g)

prev_xsol = xsol;
xsol = max(real(xsol - g), 0);
xhat = 2*xsol - prev_xsol;

rel_x = sum(power(prev_xsol(:) - xsol(:), 2));
norm_x = sum(power(xsol(:), 2));

end
