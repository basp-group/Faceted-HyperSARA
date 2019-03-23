function [y0, y, aY, sigma_noise, noise] = util_gen_input_data_noblock(im, T, W, A, input_snr)
% generates the input data


y0freq = A(im);
y0 = T * y0freq(W);


Nm = numel(y0);
normy0 = norm(y0);

Nmq = length(y0);

% add Gaussian i.i.d. noise
sigma_noise = 10^(-input_snr/20) * normy0/sqrt(Nm);

noise = (randn(Nmq, 1) + 1i*randn(Nmq, 1)) * sigma_noise/sqrt(2);
y = y0;
y = y + noise;
% normy = y./sigma_noise;
% normnoise = noise./sigma_noise;
aY = abs(y)./sigma_noise;

end





