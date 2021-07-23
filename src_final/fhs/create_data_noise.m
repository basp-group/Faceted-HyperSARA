function B = create_data_noise(y, At, G, W, Nx, Ny, No, sigma_noise, seed, operator_norm)

%! to be documented
%! create "dirty" noise matrix

rng(seed)
nChannels = numel(y);
B = zeros(Ny, Nx, nChannels);
% N = Nx*Ny;

for l = 1:nChannels
    b_l = zeros(No, 1);
    for b = 1:numel(y{l})
        noise = (randn(size(y{l}{b})) + 1i*randn(size(y{l}{b})))*sigma_noise(l)/sqrt(2);
        b_l(W{l}{b}) = b_l(W{l}{b}) + G{l}{b}' * noise;
    end
    B(:, :, l) = At(b_l);    
end

B = B/operator_norm;

end
    