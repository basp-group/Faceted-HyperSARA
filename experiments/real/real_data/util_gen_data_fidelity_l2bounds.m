function [epsilon,epsilons] = util_gen_data_fidelity_l2bounds(y_dims,Ml,param,sigma_noise)

c = length(y_dims);
epsilons = cell(numel(y_dims), 1);

if strcmp(param.type, 'sigma')
    s1 = param.sigma_ball;
    
    % estimate the global L2 bound from the chi square distribution
    epsilon = sqrt(Ml + s1*sqrt(Ml)) .* sigma_noise;
    
    % compute the data blocks L2 bounds
    for i =  1 : c
        epsilons{i} = cell(numel(y_dims{i}), 1);
        for j = 1 : numel(y_dims{i})
            m = (y_dims{i}{j});    
            epsilons{i}{j} = sqrt(m/Ml(i)) * epsilon(i);
        end
    end
end
