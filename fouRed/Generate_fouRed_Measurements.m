%% Generate Measurements

%% definition for the stoping criterion
param_l2_ball.type = 'sigma';
param_l2_ball.sigma_ball = 2;


for m = 1 : length(sigma_noise)
    [y0, y, Nm, sigma_noise_ch, noise] = util_gen_measurements_noblock(x0, G, W, A, sigma_noise(m));
    % For natural weighting
    for i = 1:length(ch)
        if normalize_data
            if numel(sigma_noise_ch{i}) == 1
                G{i} = 1./sigma_noise_ch{i} * G{i};      % Whitening G matrix (embed natural weighting in the measurement operator). In reality, this should be done by natural weighting!
                y{i} = 1./sigma_noise_ch{i} * y{i};
                noise{i} = 1./sigma_noise_ch{i} * noise{i};
            elseif numel(sigma_noise_ch{i}) == length(y{i})
                G{i} = diag(1./sigma_noise_ch{i}) * G{i};
                y{i} = diag(1./sigma_noise_ch{i}) * y{i};
                noise{i} = diag(1./sigma_noise_ch{i}) * noise{1};                    
            else
                error('Dimension of natural weights does not match')
            end
        end
    end
end

% Nm = length(ch) * numel(y_t{1}{1});
% [epsilon,epsilons] = util_gen_data_fidelity_bounds(y_t, Nm, param_l2_ball, sigma_noise);
clear y0