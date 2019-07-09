function [y0, y, Nm, sigma_noise, noise] = util_gen_measurements_block(x, G, A, input_sigma)
% generates the input data

c = size(x,3);
y_full = [];
for i = 1 : c
    Fx = A(x(:,:,i));
    for j = 1 : length(G{i})
        y0{i}{j} = G{i}{j} * Fx;
        y_full = [y_full; y0{i}{j}];
    end
end

Nm = numel(y_full(:));
normy0 = norm(y_full(:));

% add Gaussian i.i.d. noise
% sigma_noise = 10^(-input_snr/20) * normy0/sqrt(Nm);

for i = 1 : c
    for j = 1 : length(y0{i})
        sigma_noise{i}{j} = input_sigma;
        Nmq = length(y0{i}{j});
        noise{i}{j} = (randn(Nmq, 1) + 1i*randn(Nmq, 1)) * sigma_noise{i}{j};
        y{i}{j} = y0{i}{j} + noise{i}{j};
    end
end

end