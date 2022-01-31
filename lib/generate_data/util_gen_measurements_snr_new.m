function [y0, y, Ml, Nm, sigma_noise, norm_noise, G] = util_gen_measurements_snr_new(x, G, W, A, input_snr, rng_stream)
    % Generate noisy synthetic viasilities using a predefined per-channel input
    % SNR.
    %
    % Parameters
    % ----------
    % x : double[:, :, :]
    %     Input wodeband image.
    % G : cell of cell of sparse complex[:, :]
    %     Degridding matrix (per channel, per block).
    % W : cell of cell of int[:]
    %     Selection vector to extract data blocks from the full Fourier plane.
    % A : anonymous function.
    %     Function handle implementing the NUFFT.
    % input_snr : double[:]
    %     Input SNR for each spectral channel.
    % rng_stream : random number generator
    %     Random number generator.
    %
    % Returns
    % -------
    % y0 : cell of cell of complex[:]
    %     Clean (i.e., noiseless) visibilities.
    % y : cell of cell of complex[:]
    %     Noisy visibilities.
    % Ml : int[:]
    %     Number of data points in each spectral channel.
    % Nm : int
    %     Total number of visbilities (across all data blocks and channels).
    % sigma_noise : double[:]
    %     Noise level in each spectral channel.
    % norm_noise : double[:]
    %     Euclidean norm of the Gaussian noise affecting each data channel.
    % G : cell of cell of sparse complex[:, :]
    %     Degridding matrix (per channel, per block) including noise whitening
    %     contribution.
    %

    % [11/01/2022]
    c = size(x, 3);
    y0 = cell(c, 1);  % data w/o noise, not whitened
    Ml = zeros(c, 1); % number of data points per channel
    normy0 = zeros(c, 1);
    for i = 1:c
        Fx = A(x(:, :, i));
        y0{i} = cell(numel(G{i}), 1);
        for j = 1:numel(G{i})
            y0{i}{j} = G{i}{j} * Fx(W{i}{j});
            Ml(i) = numel(y0{i}{j});
            normy0(i) = normy0(i) + norm(y0{i}{j}, 'fro')^2;
        end
    end
    normy0 = sqrt(normy0);
    Nm = sum(Ml);

    % add i.i.d. Gaussian noise
    sigma_noise = (10.^(-input_snr / 20) .* normy0) ./ sqrt(Ml);

    norm_noise = zeros(c, 1);
    y = cell(c, 1);
    for i = 1:c
        y{i} = cell(numel(G{i}), 1);
        for j = 1:numel(y0{i})
            Nmq = numel(y0{i}{j});

            % previous version
            % noise = (randn(Nmq, 1) + ...
            %          1i * randn(rng_stream, Nmq, 1)) * sigma_noise(i) / sqrt(2);
            % norm_noise = norm_noise + norm(noise(:), 2)^2;
            % y{i}{j} = y0{i}{j} + noise;

            % new version (with noise whitening)
            % save non-whitened data (whitening applied in the real data code)
            noise = (randn(Nmq, 1) + ...
                     1i * randn(rng_stream, Nmq, 1)) / sqrt(2);
            y{i}{j} = y0{i}{j} + sigma_noise(i) * noise;
            norm_noise(i) = norm_noise(i) + norm(noise(:), 2)^2;
            G{i}{j} = G{i}{j} / sigma_noise(i);
        end
        norm_noise(i) = sqrt(norm_noise(i));
    end
end
