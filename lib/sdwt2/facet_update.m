% internal function
function [u_q, PsiStu_q] = facet_update(x_overlap, u_q, Iq, dims_q, I_overlap_nc_q, dims_overlap_nc_q, offset, status_q, nlevel, wavelet, N, mu, gamma_l1, rho)

soft = @(z, T) sign(z) .* max(abs(z)-T, 0);

[w, ~, ~, Ncoefs_q] = sdwt2_sara(x_overlap, Iq, dims_q, offset, status_q, nlevel, wavelet);
u_old = u_q;
w = u_q + mu * w;
u_q = w - mu * soft(w/mu, gamma_l1/mu) ;
u_q = rho*u_q + (1-rho)*u_old;
% inverse operator (for a single facet) (inverse = adjoin for spd, properly implement the adjoint operator for different boundary conditions)
PsiStu_q = isdwt2_sara(u_q-u_old, Iq, dims_q, I_overlap_nc_q, dims_overlap_nc_q, Ncoefs_q, N, nlevel, wavelet);
end