function [y, G, noise, epsilon, epsilons] = Generate_fouRed_Measurements(x0, G, A, W, input_snr, sigma_noise, usingBlocking, normalize_data)

%% Generate Measurements

%% definition for the stoping criterion
param_l2_ball.type = 'sigma';
param_l2_ball.sigma_ball = 2;


for m = 1 : length(sigma_noise)
    if usingBlocking
        [~, y, Nm, sigma_noise_ch, noise] = util_gen_measurements_block(x0, G, A, W, input_snr, sigma_noise(m)); 
    else
%     [y0, y, Nm, sigma_noise_ch, noise] = util_gen_measurements_noblock(x0, G, W, A, sigma_noise(m));
        [~, y, Nm, sigma_noise_ch, noise] = util_gen_measurements_noblock_new(x0, G, A, W, input_snr, sigma_noise(m));
    end
    % For natural weighting    
    if normalize_data
        for i = 1:length(ch)
            if usingBlocking
                for j = 1:length(G{i})
                    if numel(sigma_noise_ch{i}{j}) == 1
                        G{i}{j} = 1./sigma_noise_ch{i}{j} * G{i}{j};      % Whitening G matrix (embed natural weighting in the measurement operator). In reality, this should be done by natural weighting!
                        y{i}{j} = 1./sigma_noise_ch{i}{j} * y{i}{j};
                        noise{i}{j} = 1./sigma_noise_ch{i}{j} * noise{i}{j};
                    elseif numel(sigma_noise_ch{i}{j}) == length(y{i}{j})
                        G{i}{j} = diag(1./sigma_noise_ch{i}{j}) * G{i}{j};
                        y{i}{j} = diag(1./sigma_noise_ch{i}{j}) * y{i}{j};
                        noise{i}{j} = diag(1./sigma_noise_ch{i}{j}) * noise{i}{j};                    
                    else
                        error('Dimension of natural weights does not match')
                    end
                    [epsilon{i}{j}, epsilons{i}{j}] = util_gen_data_fidelity_bounds(y{i}{j}, Nm, param_l2_ball, 1.);
                end
            else
                if numel(sigma_noise_ch{i}) == 1
                    G{i} = 1./sigma_noise_ch{i} * G{i};      % Whitening G matrix (embed natural weighting in the measurement operator). In reality, this should be done by natural weighting!
                    y{i} = 1./sigma_noise_ch{i} * y{i};
                    noise{i} = 1./sigma_noise_ch{i} * noise{i};
                elseif numel(sigma_noise_ch{i}) == length(y{i})
                    G{i} = diag(1./sigma_noise_ch{i}) * G{i};
                    y{i} = diag(1./sigma_noise_ch{i}) * y{i};
                    noise{i} = diag(1./sigma_noise_ch{i}) * noise{i};                    
                else
                    error('Dimension of natural weights does not match')
                end
                [epsilon{i}, epsilons{i}] = util_gen_data_fidelity_bounds(y{i}, Nm, param_l2_ball, 1.);
            end
        end
        [epsilon, epsilons] = util_gen_data_fidelity_bounds(y, Nm, param_l2_ball, 1.);
    else
        [epsilon, epsilons] = util_gen_data_fidelity_bounds(y, Nm, param_l2_ball, sigma_noise(m));
    end
end

clear r u v u1 v1 uw vw
% Nm = length(ch) * numel(y_t{1}{1});
