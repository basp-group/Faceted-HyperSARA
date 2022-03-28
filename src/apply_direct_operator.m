function y = apply_direct_operator(Fx, G, Sigma)
% Apply part of the direct measurement operator (after the scaled discrete
% Fourier transform involved in the non-uniform FFT :cite:p:`Fessler2003`).
%
% Parameters
% ----------
% Fx : complex[:]
%     Scaled discrete Fourier transform of a 2d array :math:`x`.
% G : complex[:, :]
%     Gridding matrix (including noise whitening operator).
% Sigma : double[:]
%     Reduction operator.
%
% Returns
% -------
% y : complex[:]
%     Visibility vector.
%

y = Sigma .* (G * Fx);

end
