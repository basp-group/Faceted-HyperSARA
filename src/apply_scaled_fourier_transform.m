function Fx = apply_scaled_fourier_transform(x, A, No)
% Apply scaled Fourier transform to a wideband image cube :math:`x`.
%
% Parameters
% ----------
% x : double[:, :, :]
%     Input wideband image cube.
% A : anonymous function
%     anonymous function to compute the scaled Fourier transform of each slice
%     of the image
% No : int
%     Total number of pixels in the oversampled spatial Fourier plane.
%
% Returns
% -------
% Fx : double[:, :]
%     Matrix whose column contain the scaled Fourier transform of each
%     slice of ``x``.
%

Fx = zeros(No, size(x, 3));
for l = 1:size(x, 3)
    Fx(:, l) = A(x(:, :, l));
end

end
