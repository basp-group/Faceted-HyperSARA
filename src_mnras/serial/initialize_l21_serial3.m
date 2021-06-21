function [v1, weights1, s] = initialize_l21_serial3(x, Psit, extension_mode, nlevel, reweighting_alpha, sig)
%initialize_l21_serial2: initalize the dual variables related to the 
% l21-norm prior.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > x                       wideband image [M, N, L]
% > Psit                    full SARA operator @[1]                     
% > extension_mode          name of the boundary extension mode 
% > nlevel                  depth of the decomposition
%
% Output:
%
% < v1                      dual variable associated with the l21-norm 
%                           prior [s, L]
% < weights1                assciated weights for the reweigthing step 
%                           [s, L]
% < s                       number of wavelet decompostion for each channel
%                           [1]
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%

[M, N, c] = size(x);
% number of cofficients resulting from the 8 first Daubcehies wavelet
% transforms
[~, s] = n_wavelet_coefficients(2*(1:8)', [M, N], extension_mode, nlevel);
% add size from Dirac dictionary
s = s + M*N;

v1 = zeros(s, c);
for l = 1:c
    v1(:, l) = Psit(x(:, :, l));
end
% compute part of the row-wise l2 norm
d1_ = sum(v1.^2, 2); 
% compute sum across all workers
d1 = gplus(d1_);
d1 = sqrt(d1);
upsilon = sig*reweighting_alpha;
weights1 = upsilon ./ (upsilon + d1);

v1 = zeros(s, c);

end
        