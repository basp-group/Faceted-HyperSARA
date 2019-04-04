%% Generate Measurements

%% definition for the stoping criterion
param_l2_ball.type = 'sigma';
param_l2_ball.sigma_ball = 2;


%%
q = 1;
for m = 1 : length(sigma_noise)
    for n = 1 : num_tests
        
        [y0, y, Nm] = util_gen_measurements(x0, G, W, A, sigma_noise(m),seed);
        
        [epsilon,epsilons] = util_gen_data_fidelity_bounds(y, Nm, param_l2_ball, sigma_noise(m));
        
        y0_t{q} = y0;
        y_t{q} = y;
        epsilon_t{q} = epsilon;
        epsilons_t{q} = epsilons;
        
        q = q + 1;
        
    end
end