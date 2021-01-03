function [v0, v1, weights0, weights1] = initialize_dual_variables_prior_overlap(Ncoefs, dims, dims_overlap_ref, dirac_present, c, nlevel)

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
