function [y0, y, Nm] = util_gen_measurements_sigma(x, G, W, A, sigma_noise, seed)
% [extended_summary]
% 
% Parameters
% ----------
% x : [type]
%     [description]
% G : [type]
%     [description]
% W : [type]
%     [description]
% A : [type]
%     [description]
% sigma_noise : [type]
%     [description]
% seed : [type]
%     [description]
% 
% Returns
% -------
% [type]
%     [description]

    % generates the input data
    rng(seed);

    c = size(x, 3);
    y_full = [];
    for i = 1:c
        Fx = A(x(:, :, i));
        for j = 1:length(G{i})
            y0{i}{j} = G{i}{j} * Fx(W{i}{j});
            y_full = [y_full; y0{i}{j}];
        end
    end

    Nm = numel(y_full(:));
    normy0 = norm(y_full(:));
    % sigma_noise = 0.1

    % add Gaussian i.i.d. noise
    % sigma_noise = 10^(-input_snr/20) * normy0/sqrt(Nm);

    norm_noise = 0;
    for i = 1:c
        for j = 1:length(y0{i})
            Nmq = length(y0{i}{j});
            % noise = (randn(Nmq, 1) + 1i*randn(Nmq, 1)) * sigma_noise;
            noise = (randn(Nmq, 1) + 1i * randn(Nmq, 1)) * sigma_noise / sqrt(2);
            norm_noise = norm_noise + norm(noise(:), 2)^2;
            y{i}{j} = y0{i}{j} + noise;
        end
    end
    norm_noise = sqrt(norm_noise);

end
