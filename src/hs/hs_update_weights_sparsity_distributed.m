function weights = hs_update_weights_sparsity_distributed(x, Psit, ...
    weights, reweight_alpha, sig)
% Update the weigths for the reweighting of the
% l21-norm prior.
%
% Parameters
% ----------
% x : array (3d)
%     Wideband image [M, N, L].
% Psit : lambda functions
%     Full SARA operator @[1].
% weights : array (2d)
%     Weights associated with the reweigthing step [s, L].
% reweight_alpha : double
%     Reweighting parameter [1].
% sig : double
%     Estimate of the noise level in the SARA domain.
%
% Returns
% -------
% weights : array (2d)
%     Updated weights associated with the reweigthing step [s, L].
%

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%

w = zeros(size(weights, 1), size(x, 3));
for l = 1 : size(x, 3)
    w(:, l) = Psit(x(:, :, l));
end
d_ = sum(abs((w)).^2,2);
d = gplus(d_);
d = sqrt(d);
upsilon = sig*reweight_alpha;
weights = upsilon ./ (upsilon + d);

end
