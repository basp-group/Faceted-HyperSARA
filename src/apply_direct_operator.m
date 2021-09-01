function y = apply_direct_operator(Fx, G, Sigma)
% Apply part of the direct measurement operator (after the scaled discrete 
% Fourier transform involved in the non-uniform FFT :cite:p:`Fessler2003`).
%
% Parameters
% ----------
% Fx : array (complex)
%     Scaled discrete Fourier transform of a 2d array :math:`x`.
% G : array (real)
%     Gridding matrix (including noise whitening operator).
% Sigma : array (real)
%     Reduction operator.
%
% Returns
% -------
% y : array (complex)
%     Visibility vector.
%

y = Sigma .* (G * Fx);

end