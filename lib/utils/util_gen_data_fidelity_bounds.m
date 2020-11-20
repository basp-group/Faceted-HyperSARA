function [epsilon,epsilons] = util_gen_data_fidelity_bounds(y,Nm,param,sigma_noise)

c = length(y);

if param.type == 'sigma'
    s1 = param.sigma_ball;
    
    % estimate the global L2 bound from the chi square distribution
    epsilon = sqrt(Nm + s1*sqrt(Nm)) * sigma_noise
    
    % compute the data blocks L2 bounds
    for i =  1 : c
        for j = 1 : length(y{i})
        m = length(y{i}{j});    
        epsilons{i}{j} = sqrt(m/Nm) * epsilon;
        end
    end
end









