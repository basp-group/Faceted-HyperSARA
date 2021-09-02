function [global_norm_res_c, norm_epsilon_c, global_norm_res_a, norm_epsilon_a] = ...
    sanity_check_rd(epsilon, norm_res)
% Update the data fidelity term in the preconditioned primal-dual algorithm.
%
% Update the data fidelity terms owned by each worked involved in the group
% of data nodes (with preconditioning -> projection onto an ellipsoid).
%
% Args:
%     norm_res (cell): norm of the residual {L}{nblocks}[1].
%     epsilon (cell): l2-ball constraints {L}{nblocks}[1].
%
% Returns:
%     global_norm_res (double): square global norm of the residual.
%     norm_epsilon (double): square global norm of epsilon.

% -------------------------------------------------------------------------%
%%
% Reference:
%
% [1] Dabbech, A. and Onose, A. and Abdulaziz, A. and Perley, R. A. and
% Smirnov, O. M. and Wiaux, Y.Cygnus A Super-Resolved via Convex
% Optimization from VLA Data, Mon. Not. R. Astron. Soc., 476(3),
% pp. 2853-2866
% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [...]
% -------------------------------------------------------------------------%
%%

nChannels = numel(norm_res);
norm_epsilon_c = 0;
global_norm_res_c = 0;
norm_epsilon_a = 0;
global_norm_res_a = 0;
for i = 1:nChannels
    for j = 1:numel(norm_res{i})
        if j == 1  %%%% C block
            global_norm_res_c = global_norm_res_c + norm_res{i}{j}^2;
            norm_epsilon_c = norm_epsilon_c + power(epsilon{i}{j}, 2);
        else %%%% A block
            global_norm_res_a = global_norm_res_a + norm_res{i}{j}^2;
            norm_epsilon_a = norm_epsilon_a + power(epsilon{i}{j}, 2);
        end
    end
end

end
