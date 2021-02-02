function [epsilon, t_block] = update_epsilon(epsilon, t, t_block, rel_fval,...
    norm_res, adapt_eps_tol_in, adapt_eps_tol_out, adapt_eps_steps, ...
    adapt_eps_change_percentage, l2_upper_bound)
% Adaptive update of the l2-ball constraint.
%
% Update of the l2-ball constraints of radius epsilon (adaptive 
% :math:`\varepsilon`), based on the scheme introduced in :cite:`Dabbech2018`).
%
% Args:
%     epsilon (cell): l2-ball radii {L}{nblocks}.
%     t (int): current iteration index.
%     t_block (cell): last iteration at which each block has
%                     been updated {L}{nblocks}.
%     rel_fval (double): relative variation of the solution.
%     norm_res (cell): norm of the residual image {L}{nblocks}.
%     adapt_eps_tol_in (double): tolerance for the norm of the residual 
%                                (below the current value of epsilon,
%                                within the ball).
%     adapt_eps_tol_out (double): tolerance for the norm of the residual 
%                                 (above current value of epsilon, out of
%                                 the ball).
%     adapt_eps_steps (double): number of steps between two consecutive
%                               updates.
%     adapt_eps_rel_obj (double): relative variation of the objective
%                                 function
%     adapt_eps_change_percentage (double): update factor for epsilon.
%
% Returns:
%     epsilon (cell): updated l2-ball radii {L}{nblocks}.
%     t_block (cell): last iteration at which each block has been updated 
%                     {L}{nblocks}.
%

%-------------------------------------------------------------------------%
%%
% Reference:
%
% [1] Dabbech, A. and Onose, A. and Abdulaziz, A. and Perley, R. A. and 
% Smirnov, O. M. and Wiaux, Y. Cygnus A Super-Resolved via Convex 
% Optimization from VLA Data, Mon. Not. R. Astron. Soc., 476(3), 
% pp. 2853-2866
%
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%
nChannels = length(epsilon);

% if rel_fval < adapt_eps_rel_obj % can be done out of the function
for i = 1 : nChannels
    for  j = 1 : length(epsilon{i})
        if t > t_block{i}{j} + adapt_eps_steps %&& rel_fval < adapt_eps_rel_obj
            if  norm_res{i}{j} < adapt_eps_tol_in * epsilon{i}{j}
                epsilon{i}{j} = adapt_eps_change_percentage*norm_res{i}{j} + (1 - adapt_eps_change_percentage)*epsilon{i}{j};
                t_block{i}{j} = t;
                fprintf('Updated  epsilon DOWN: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j},norm_res{i}{j},j,i);
            end

            if norm_res{i}{j} > adapt_eps_tol_out * epsilon{i}{j}
                target_eps = adapt_eps_change_percentage*norm_res{i}{j} + (1 - adapt_eps_change_percentage)*epsilon{i}{j};
                if target_eps > l2_upper_bound{i}{j}
                    epsilon{i}{j} = l2_upper_bound{i}{j};
                else
                   epsilon{i}{j} = target_eps;
                end
                t_block{i}{j} = t;
                fprintf('Updated  epsilon UP: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j},norm_res{i}{j},j,i);
            end
        end
    end
end
% end

end
