function weights = hs_update_weights_lowrankness_serial(x, ...
    reweight_alpha, sig_bar)
% Update the weigths for the reweighting of the nuclear norm prior.
%
% Parameters
% ----------
% x : array (3d)
%     Wideband image [M, N, L].
% reweight_alpha : double
%     Reweighting parameter [1].
% sig_bar : double
%     Noise level (singular value space).
%
% Returns
% -------
% weights : array (2d)
%     Weights associated with the reweigthing step [s, L].
%
             
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%

[M,  N, c] = size(x);
[~, D, ~] = svd(reshape(x, [M*N, c]),'econ');
d = abs(diag(D));
upsilon_bar = sig_bar*reweight_alpha;
weights = upsilon_bar ./ (upsilon_bar + d);

end
