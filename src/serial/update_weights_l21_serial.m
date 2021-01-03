function weights = update_weights_l21_serial(x, Psit, weights, reweight_alpha, sig)
%update_weights_l21_serial: update the weigths for the reweighting of the
% l21-norm prior.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x                       wideband image [M, N, L]
% > Psit                    full SARA operator @[1]                     
% > weights                 weights associated with the reweigthing step 
%                           [s, L]
% > reweight_alpha          reweighting parameter [1]
%
% Output:
%
% < weights                 weights associated with the reweigthing step 
%                           [s, L]
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
d = sqrt(sum(abs((w)).^2,2));
upsilon = sig*reweight_alpha;
weights = upsilon ./ (upsilon + d);

end
