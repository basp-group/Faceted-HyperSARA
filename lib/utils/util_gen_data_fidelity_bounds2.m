function [epsilon,epsilons] = util_gen_data_fidelity_bounds2(y,Ml,param,sigma_noise)

c = length(y);
epsilons = cell(numel(y), 1);

if strcmp(param.type, 'sigma')
    s1 = param.sigma_ball;
    
    % estimate the global L2 bound from the chi square distribution
    epsilon = sqrt(Ml + s1*sqrt(Ml)) .* sigma_noise;
    
    % compute the data blocks L2 bounds
    for i =  1 : c
        epsilons{i} = cell(numel(y{i}), 1);
        for j = 1 : numel(y{i})
            m = numel(y{i}{j});    
            epsilons{i}{j} = sqrt(m/Ml(i)) * epsilon(i);
        end
    end
end
