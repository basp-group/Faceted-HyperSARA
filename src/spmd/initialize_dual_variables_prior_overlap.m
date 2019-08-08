function [v0, v1, weights0, weights1] = ...
    initialize_dual_variables_prior_overlap(Ncoefs, dims, dims_overlap_ref,...
                                            dirac_present, c, nlevel)
%initialize_dual_variables_prior_overlap: initialize all the dual variables
% for a given facet (nuclear and l21 norms).
%-------------------------------------------------------------------------%
%%
% Input:
%
% > Ncoefs              number of wavelet coefficients at each scale
% > dims                size of a non-overlapping facet [1, 2]
% > dims_overlap_ref    dimension of a facet (with overlap) [1, 2]
% > dirac_present       flag indicating whether the Dirac dictionary is
%                       used in the sparsifying dictionary
% > c                   number of spectral channels [1]
% > nlevel              number of wavelet decompositions [1]
%
% Output:
%
% < v0                  dual variable associated with the nuclear norm
% < v1                  dual variable associated with the l21-norm
% < weights0            weigths associated with the nuclear norm
% < weights1            weigths ssociated with the l21-norm
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
%-------------------------------------------------------------------------%
%%

% compute size of the dual variables and the weigths
p = prod(Ncoefs, 2);
% if dirac_present
%     sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) + prod(dims);
% else
%     sz = 3*sum(p) - 2*sum(p(nlevel+1:nlevel+1:end));
% end
sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) - 2*p(end); % number of coeffs with the Dirac basis

v0 = zeros(prod(dims_overlap_ref), c);
weights0 = ones(min(prod(dims_overlap_ref), c), 1);
v1 = zeros(sz, c);
weights1 = ones(sz, 1);

end
