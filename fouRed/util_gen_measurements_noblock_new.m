function [y0, y, Nm, sigma_noise, noise] = util_gen_measurements_noblock_new(x, G, A, input_sigma)

c = size(x,3);
y_full = [];
for i = 1 : c
    Fx = A(x(:,:,i));
    y0{i} = G{i} * Fx;
    y_full = [y_full; y0{i}];
end

Nm = numel(y_full(:));
normy0 = norm(y_full(:));

% add Gaussian i.i.d. noise
sigma_noise = cell(1, c);

for i = 1:c
%     sigma_noise{i} = 10^(-input_snr/20) * normy0/sqrt(Nm);
    sigma_noise{i} = input_sigma;
end

for i = 1 : c
    Nmq = length(y0{i});
    noise{i} = (randn(Nmq, 1) + 1i*randn(Nmq, 1)) * sigma_noise{i};
%     noise{i} = (randn(Nmq, 1) + 1i*randn(Nmq, 1)) * sigma_noise{i}/sqrt(2);
    y{i} = y0{i} + noise{i};
end

end