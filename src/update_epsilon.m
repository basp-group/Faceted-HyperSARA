function [epsilon, t_block] = update_epsilon(epsilon, t, t_block, rel_fval, norm_res, ...
    adapt_eps_tol_in, adapt_eps_tol_out, adapt_eps_steps, adapt_eps_rel_obj, adapt_eps_change_percentage)

ci = length(epsilon);

% to be rewritten in a simpler manner
for i = 1 : ci
    for  j = 1 : length(epsilon{i})
        if  norm_res{i}{j} < adapt_eps_tol_in * epsilon{i}{j}
            if t > t_block{i}{j} + adapt_eps_steps && rel_fval(t) < adapt_eps_rel_obj
                epsilon{i}{j} = norm_res{i}{j} + (-norm_res{i}{j} + epsilon{i}{j}) * (1 - adapt_eps_change_percentage);
                t_block{i}{j} = t;
                %count_eps_update_down = count_eps_update_down + 1;
                %fprintf('Updated  epsilon DOWN: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j},norm_res{i}{j},j,i);
            end
        end
        
        if  norm_res{i}{j} > adapt_eps_tol_out * epsilon{i}{j}
            if t > t_block{i}{j} + adapt_eps_steps && rel_fval(t) < adapt_eps_rel_obj
                epsilon{i}{j} = epsilon{i}{j} + (norm_res{i}{j} - epsilon{i}{j}) * adapt_eps_change_percentage;
                t_block{i}{j} = t;
                %count_eps_update_up = count_eps_update_up + 1;
                %fprintf('Updated  epsilon UP: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j},norm_res{i}{j},j,i);
            end
        end
    end
end

end