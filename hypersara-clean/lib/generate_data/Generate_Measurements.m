%% Generate Measurements

%% definition for the stoping criterion
param_l2_ball.type = 'sigma';
param_l2_ball.sigma_ball = 2;


%%
[y0, y, Nm] = util_gen_measurements(x0, G, W, A, sigma_noise,seed);

[epsilon,epsilons] = util_gen_data_fidelity_bounds(y, Nm, param_l2_ball, sigma_noise);