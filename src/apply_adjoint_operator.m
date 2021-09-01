function Fx = apply_adjoint_operator(y, G, Sigma)
% Apply part of the adjoint measurment operator (before the scaled 
% inverse discrete Fourier transform involved in the non-uniform FFT 
% :cite:p:`Fessler2003`).
%
% Parameters
% ----------
% y : array (complex)
%     Input data array (data-block).
% G : sparse array
%     Gridding matrix (including noise whitening operator).
% Sigma : array (real)
%     Reduction operator.
%
% Returns
% -------
% Fx : array (complex)
%     Gridded visibilities.
%
    
Fx = G' * (Sigma .* y);

end