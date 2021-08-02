function Fx = apply_scaled_fourier(x, A, No)
% Apply scaled Fourier transform to a wideband image cube.
%
% Parameters
% ----------
% x : array (3d)
%     Input wideband image cube.
% A : lambda function
%     Lambda function to compute the scaled Fourier transform of each slice
%     of the image
% No : int
%     Total number of pixels in the oversampled spatial Fourier plane.
%
% Returns
% -------
% Fx : array (2d)
%     Matrix whose column contain the scaled Fourier transform of each 
%     slice of ``x``.
%

Fx = zeros(No, size(x, 3));
for l = 1:size(x, 3)
    Fx(:,l) = A(x(:, :, l));
end

end