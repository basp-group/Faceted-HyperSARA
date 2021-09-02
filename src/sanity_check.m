function [global_norm_res, norm_epsilon] = sanity_check(epsilon, norm_res)
% Update the data fidelity term in the preconditioned primal-dual algorithm.
%
% Update the data fidelity terms owned by each worked involved in the group
% of data nodes (with preconditioning -> projection onto an ellipsoid).
%
% Parameters
% ----------
% epsilon : cell
%     Norm of the residual {L}{nblocks}[1].
% norm_res : cell
%     :math:`\ell_2`-ball constraints {L}{nblocks}[1].
%
% Returns
% -------
% global_norm_res : double
%     Square global norm of the residual.
% norm_epsilon : double
%     Square global norm of epsilon.
%

%%
nChannels = numel(norm_res);
norm_epsilon = 0;
global_norm_res = 0;
for i = 1:nChannels
    for j = 1:numel(norm_res{i})
        global_norm_res = global_norm_res + norm_res{i}{j}^2;
        norm_epsilon = norm_epsilon + power(epsilon{i}{j}, 2);
    end
end

end
