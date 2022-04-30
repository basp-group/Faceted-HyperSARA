function Fx = apply_adjoint_operator(y, G, Lambda)
% Apply part of the adjoint measurment operator (before the scaled
% inverse discrete Fourier transform involved in the non-uniform FFT
% :cite:p:`Fessler2003`).
%
% Parameters
% ----------
% y : complex[:]
%     Input data array (data-block).
% G : sparse array
%     Gridding matrix (including noise whitening operator).
% Lambda : double[:]
%     Reduction operator.
%
% Returns
% -------
% Fx : complex[:]
%     Gridded visibilities.
%

Fx = G' * (Lambda .* y);

end
