%% Generate Measurements

%% definition for the stoping criterion
param_l2_ball.type = 'sigma';
param_l2_ball.sigma_ball = 2;


for m = 1 : length(sigma_noise)
    if usingBlocking
        [~, y, Nm, sigma_noise_ch, noise] = util_gen_measurements_block(x0, Gw, A, sigma_noise(m));                
    else
%     [y0, y, Nm, sigma_noise_ch, noise] = util_gen_measurements_noblock(x0, G, W, A, sigma_noise(m));
        [~, y, Nm, sigma_noise_ch, noise] = util_gen_measurements_noblock_new(x0, Gw, A, sigma_noise(m));
    end
    % For natural weighting    
    if normalize_data
        for i = 1:length(ch)
            if usingBlocking
                for j = 1:length(Gw{i})
                    if numel(sigma_noise_ch{i}{j}) == 1
                        Gw{i}{j} = 1./sigma_noise_ch{i}{j} * Gw{i}{j};      % Whitening G matrix (embed natural weighting in the measurement operator). In reality, this should be done by natural weighting!
                        y{i}{j} = 1./sigma_noise_ch{i}{j} * y{i}{j};
                        noise{i}{j} = 1./sigma_noise_ch{i}{j} * noise{i}{j};
                    elseif numel(sigma_noise_ch{i}{j}) == length(y{i}{j})
                        Gw{i}{j} = diag(1./sigma_noise_ch{i}{j}) * Gw{i}{j};
                        y{i}{j} = diag(1./sigma_noise_ch{i}{j}) * y{i}{j};
                        noise{i}{j} = diag(1./sigma_noise_ch{i}{j}) * noise{i}{j};                    
                    else
                        error('Dimension of natural weights does not match')
                    end
                    
                end
            else
                if numel(sigma_noise_ch{i}) == 1
                    Gw{i} = 1./sigma_noise_ch{i} * Gw{i};      % Whitening G matrix (embed natural weighting in the measurement operator). In reality, this should be done by natural weighting!
                    y{i} = 1./sigma_noise_ch{i} * y{i};
                    noise{i} = 1./sigma_noise_ch{i} * noise{i};
                elseif numel(sigma_noise_ch{i}) == length(y{i})
                    Gw{i} = diag(1./sigma_noise_ch{i}) * Gw{i};
                    y{i} = diag(1./sigma_noise_ch{i}) * y{i};
                    noise{i} = diag(1./sigma_noise_ch{i}) * noise{i};                    
                else
                    error('Dimension of natural weights does not match')
                end
            end
        end
    end
end

clear r u v u1 v1 uw vw
% Nm = length(ch) * numel(y_t{1}{1});
% [epsilon,epsilons] = util_gen_data_fidelity_bounds(y_t, Nm, param_l2_ball, sigma_noise);