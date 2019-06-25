function [epsilon, t_block] = update_epsilon(epsilon, t, t_block, rel_fval,...
    norm_res, adapt_eps_tol_in, adapt_eps_tol_out, adapt_eps_steps, ...
    adapt_eps_rel_obj, adapt_eps_change_percentage)
%update_epsilon: adaptive epsilon scheme (i.e., update of the l2-ball radii
% epsilon, based on the technique introduced in [1]).
%-------------------------------------------------------------------------%
%%
% Reference:
%
% [1] Dabbech, A. and Onose, A. and Abdulaziz, A. and Perley, R. A. and 
% Smirnov, O. M. and Wiaux, Y.Cygnus A Super-Resolved via Convex 
% Optimization from VLA Data, Mon. Not. R. Astron. Soc., 476(3), 
% pp. 2853-2866
%-------------------------------------------------------------------------%
%%
% Input:
%
% > epsilon                         l2-ball radii {L}{nblocks}
% > t                               current iteration index [1]
% > t_block                         last iteration at which each block has
%                                   been updated {L}{nblocks}
% > rel_fval                        relative variation of the solution [1]
% > norm_res                        norm of the residual image {L}{nblocks}
% > adapt_eps_tol_in                ...
% > adapt_eps_tol_out               ...
% > adapt_eps_steps                 number of steps between two consecutive
%                                   updates [1]
% > adapt_eps_rel_obj               relative variation of the objective
%                                   function [1]
% > adapt_eps_change_percentage     update factor for epsilon [1]
%
% Output:
%
% < epsilon                         updated l2-ball radii {L}{nblocks}
% < t_block                         last iteration at which each block has
%                                   been updated {L}{nblocks}
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [../../2019]
%-------------------------------------------------------------------------%
%%

nChannels = length(epsilon);

for i = 1 : nChannels
    for  j = 1 : length(epsilon{i})
        if  norm_res{i}{j} < adapt_eps_tol_in * epsilon{i}{j}
            if t > t_block{i}{j} + adapt_eps_steps && rel_fval < adapt_eps_rel_obj
                epsilon{i}{j} = norm_res{i}{j} + (-norm_res{i}{j} + epsilon{i}{j}) * (1 - adapt_eps_change_percentage);
                t_block{i}{j} = t;
                fprintf('Updated  epsilon DOWN: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j},norm_res{i}{j},j,i);
            end
        end
        
        if  norm_res{i}{j} > adapt_eps_tol_out * epsilon{i}{j}
            if t > t_block{i}{j} + adapt_eps_steps && rel_fval < adapt_eps_rel_obj
                epsilon{i}{j} = epsilon{i}{j} + (norm_res{i}{j} - epsilon{i}{j}) * adapt_eps_change_percentage;
                t_block{i}{j} = t;
                fprintf('Updated  epsilon UP: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j},norm_res{i}{j},j,i);
            end
        end
    end
end

end
