% internal function
function [PsiStu_q, u_q] = facet_update_par(x_overlap, u_q, Iq, dims_q, I_overlap_nc_q, dims_overlap_nc_q, offsetp, status_q, nlevelp, waveletp, Np, mup, gamma_l1p, rhop)

soft = @(z, T) sign(z) .* max(abs(z)-T, 0);

[w, ~,  ~, Ncoefs_q] = sdwt2_sara(x_overlap, Iq, dims_q, offsetp.Value, status_q, nlevelp.Value, waveletp.Value);
u_old = u_q;
w = u_q + mup.Value * w;
u_q = w - mup.Value * soft(w/mup.Value, gamma_l1p.Value/mup.Value) ;
u_q = rhop.Value*u_q + (1-rhop.Value)*u_old;
% inverse operator (for a single facet) (inverse = adjoin for zpd, properly implement the adjoint operator for different boundary conditions)
PsiStu_q = isdwt2_sara(u_q-u_old, Iq, dims_q, I_overlap_nc_q, dims_overlap_nc_q, Ncoefs_q, Np.Value, nlevelp.Value, waveletp.Value);
end