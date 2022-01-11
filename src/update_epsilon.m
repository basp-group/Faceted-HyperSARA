function [epsilon, t_block] = update_epsilon(epsilon, t, t_block, ...
                                            norm_res, adapt_eps_tol_in, ...
                                            adapt_eps_tol_out, ...
                                            adapt_eps_steps, ...
                                            adapt_eps_change_percentage)
% Adaptive update of the :math:`\ell_2`-ball constraint.
%
% Update of the :math:`\ell_2`-ball constraints of radius ``epsilon``
% (adaptive :math:`\varepsilon`), based on the scheme introduced in
% :cite:p:`Dabbech2018`).
%
% Parameters
% ----------
% epsilon : cell
%     Radius of the :math:`\ell_2`-ball constraint {L}{nblocks}.
% t : int
%     Current iteration index.
% t_block : cell
%     Last iteration at which each block of ``epsilon`` constraint has been
%     updated {L}{nblocks}.
% norm_res : cell
%     Norm of the residual image {L}{nblocks}.
% adapt_eps_tol_in : double
%     Tolerance for the norm of the residual (below the current value of
%     ``epsilon``, within the ball).
% adapt_eps_tol_out : double
%     Tolerance for the norm of the residual (above current value of
%     ``epsilon``, out of the ball).
% adapt_eps_steps : int
%     Number of iteration between two consecutive updates.
% adapt_eps_change_percentage : double
%     Update factor for ``epsilon``.

% Returns
% -------
% epsilon : cell
%     Updated :math:`\ell_2`-ball radii {L}{nblocks}.
% t_block : cell
%     Last iteration at which each block has been updated {L}{nblocks}.
%

% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
% TODO: update name of variable
% -------------------------------------------------------------------------%

%%
n_channels = length(epsilon);

for i = 1:n_channels
    for  j = 1:length(epsilon{i})
        if t > t_block{i}{j} + adapt_eps_steps
            if  norm_res{i}{j} < adapt_eps_tol_in * epsilon{i}{j}
                epsilon{i}{j} = adapt_eps_change_percentage * norm_res{i}{j} + ...
                    (1 - adapt_eps_change_percentage) * epsilon{i}{j};
                t_block{i}{j} = t;
                fprintf('Updated  epsilon DOWN: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j}, norm_res{i}{j}, j, i);
            end

            if norm_res{i}{j} > adapt_eps_tol_out * epsilon{i}{j}
                target_eps = adapt_eps_change_percentage * norm_res{i}{j} + ...
                    (1 - adapt_eps_change_percentage) * epsilon{i}{j};
                epsilon{i}{j} = target_eps;
                t_block{i}{j} = t;
                fprintf('Updated  epsilon UP: %e\t, residual: %e\t, Block: %i, Band: %i\n', epsilon{i}{j}, norm_res{i}{j}, j, i);
            end
        end
    end
end

end
