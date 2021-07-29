function Fx = apply_scaled_fourier(x, A, No)

Fx = zeros(No, size(x, 3));
for l = 1:size(x, 3)
    Fx(:,l) = A(x(:, :, l));
end

end