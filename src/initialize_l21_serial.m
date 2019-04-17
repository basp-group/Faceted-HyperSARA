function [v1, weights1, s] = initialize_l21_serial(x, Psit, extension_mode, nlevel)

[M, N, c] = size(x);
[~, s] = n_wavelet_coefficients(2*(1:8)', [M, N], extension_mode, nlevel);
s = s + M*N; % Dirac dictionary

v1 = zeros(s, c);
weights1 = ones(s, 1);

for l = 1:c
    v1(:, l) = Psit(x(:, :, l));
end

end
