function [v1, weights1, s] = initialize_l21_serial(x, Psit, extension_mode, nlevel)
%run_par_l21: update the dual variable related to the facet l21-norm
% prior.
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
% [../../2019]
%-------------------------------------------------------------------------%
%%

[M, N, c] = size(x);
[~, s] = n_wavelet_coefficients(2*(1:8)', [M, N], extension_mode, nlevel);
s = s + M*N; % Dirac dictionary

v1 = zeros(s, c);
weights1 = ones(s, 1);

for l = 1:c
    v1(:, l) = Psit(x(:, :, l));
end

end
