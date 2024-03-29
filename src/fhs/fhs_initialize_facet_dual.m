function [v0, v1, weights0, weights1] = fhs_initialize_facet_dual(Ncoefs, ...
    dims_o, n_channels, nlevel)
% Initialize dual variables (constant overlap).
%
% Initialize all the dual variables for a given facet (nuclear and
% :math:`\ell_{2, 1}` norms).
%
% Parameters
% ----------
% Ncoefs : int[:, :]
%     Number of wavelet coefficients at each scale.
% dims_o :int[:]
%     Dimension of a facet (with overlap) ``[1, 2]``.
% n_channels : int
%     Number of spectral channels.
% nlevel : int
%     Depth of the wavelet decompositions.
%
% Returns
% -------
% v0 : double[:, :]
%     Dual variable associated with the nuclear norm.
% v1 : double[:, :]
%     Dual variable associated with the :math:`\ell_{2,1}` norm.
% weights0 : double[:]
%     Weigths associated with the nuclear norm.
% weights1 : double[:]
%     Weigths ssociated with the :math:`\ell_{2,1}` norm.
%

% -------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [08/08/2019]
% -------------------------------------------------------------------------%
%%

% compute size of the dual variables and the weigths
p = prod(Ncoefs, 2);
% if dirac_present
%     sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) + prod(dims);
% else
%     sz = 3*sum(p) - 2*sum(p(nlevel+1:nlevel+1:end));
% end

% number of SARA coeffs with the Dirac basis
sz = 3 * sum(p(1:end)) - 2 * sum(p(nlevel + 1:nlevel + 1:end)) - 2 * p(end);

v0 = zeros(prod(dims_o), n_channels);
weights0 = ones(min(prod(dims_o), n_channels), 1);
v1 = zeros(sz, n_channels);
weights1 = ones(sz, 1);

end
