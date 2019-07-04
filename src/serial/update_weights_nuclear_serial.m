function weights = update_weights_nuclear_serial(x, reweight_alpha)
%update_weights_nuclear_serial: update the weigths for the reweighting of 
% the nuclear norm prior.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x                       wideband image [M, N, L]
% > reweight_alpha          reweighting parameter [1]
%
% Output:
%
% < weights                 weights associated with the reweigthing step 
%                           [s, L]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% [../../2019]
%-------------------------------------------------------------------------%
%%

[M,  N, c] = size(x);
[~, D, ~] = svd(reshape(x, [M*N, c]),'econ');
d = abs(diag(D));
weights = reweight_alpha ./ (reweight_alpha + d);

end
