function [epsilon, t_block] = update_epsilon(epsilon, t, t_block, ...
    norm_res, adapt_eps_tol_in, adapt_eps_tol_out, adapt_eps_steps, ...
    adapt_eps_change_percentage) % l2_upper_bound
% Adaptive update of the :math:`\ell_2`-ball constraint.
%
% Update of the :math:`\ell_2`-ball constraints of radius ``epsilon`` 
% (adaptive :math:`\varepsilon`), based on the scheme introduced in
% :cite:p:`Dabbech2018`).
%
% Parameters
% ----------
% epsilon : cell
%     :math:`\ell_2`-ball constraint {L}{nblocks}.
% t : int
%     current iteration index.
% t_block : cell
%     (cell): last iteration at which each block has been updated {L}{nblocks}.
% norm_res : cell
%     norm of the residual image {L}{nblocks}.
% adapt_eps_tol_in : double
%     tolerance for the norm of the residual (below the current value of 
%     ``epsilon``, within the ball).
% adapt_eps_tol_out : double
%     tolerance for the norm of the residual (above current value of 
%     ``epsilon``, out of the ball).
% adapt_eps_steps : int
%     umber of steps between two consecutive updates.
% adapt_eps_change_percentage : double
%     update factor for ``epsilon``.

% Returns
% -------
% epsilon : cell
%     updated :math:`\ell_2`-ball radii {L}{nblocks}.
% t_block : cell
%     last iteration at which each block has been updated {L}{nblocks}.
%

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
% TODO: update name of variable
%-------------------------------------------------------------------------%
%%
nChannels = length(epsilon);

for i = 1 : nChannels
    for  j = 1 : length(epsilon{i})
        if t > t_block{i}{j} + adapt_eps_steps
            if  norm_res{i}{j} < adapt_eps_tol_in * epsilon{i}{j}
                epsilon{i}{j} = adapt_eps_change_percentage*norm_res{i}{j} + (1 - adapt_eps_change_percentage)*epsilon{i}{j};
                t_block{i}{j} = t;
                fprintf('Updated  epsilon DOWN: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j},norm_res{i}{j},j,i);
            end

            if norm_res{i}{j} > adapt_eps_tol_out * epsilon{i}{j}
                target_eps = adapt_eps_change_percentage*norm_res{i}{j} + (1 - adapt_eps_change_percentage)*epsilon{i}{j};
                epsilon{i}{j} = target_eps;
                t_block{i}{j} = t;
                fprintf('Updated  epsilon UP: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j},norm_res{i}{j},j,i);
            end
        end
    end
end

end
