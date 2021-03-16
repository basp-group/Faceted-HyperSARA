function [v0, v1, weights0, weights1] = initialize_dual_and_weights2(x_overlap, ...
    I, offset, status, nlevel, wavelet, Ncoefs, dims_overlap_ref, ...
    offsetL, offsetR, reweight_alpha, crop_l21, crop_nuclear, w)
% Initialize dual variables (constant overlap). 
% ![DOCUMENTATION TO BE UPDATED]
%
% Initialize all the dual variables for a given facet (nuclear and l21 
% norms).
%
% Args:
%     Ncoefs (array): number of wavelet coefficients at each scale.
%     dims_o (array): dimension of a facet (with overlap) [1, 2].
%     c (int): number of spectral channels.
%     nlevel (int): depth of decomposition.
%
% Returns:
%     v0 (array): dual variable associated with the nuclear norm.
%     v1 (array): dual variable associated with the l21-norm.
%     weights0 (array): weigths associated with the nuclear norm.
%     weights1 (array): weigths ssociated with the l21-norm.

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [../../....]
%-------------------------------------------------------------------------%
%%

%! dual variables and weights initializaed from the current value of 
%! the primal variable

% compute size of the dual variables and the weigths
p = prod(Ncoefs, 2);
% if dirac_present
%     sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) + prod(dims);
% else
%     sz = 3*sum(p) - 2*sum(p(nlevel+1:nlevel+1:end));
% end
sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) - 2*p(end); % number of coeffs with the Dirac basis

% nuclear norm
% weights0 = ones(min(prod(dims_o), c), 1);
% v0 = zeros(prod(dims_o), c);
v0 = w.*x_overlap(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :);
sol = reshape(v0, [numel(v0)/size(v0, 3), size(x_overlap, 3)]);
[~,S00,~] = svd(sol,'econ');
d0 = abs(diag(S00));
% upsilon_bar = sig_bar*reweight_alpha;
% weights0 = upsilon_bar ./ (upsilon_bar + d0);
weights0 = reweight_alpha ./ (reweight_alpha + d0);

% l21 norm
zerosNum = dims_overlap_ref + offsetL + offsetR; % offset for the zero-padding
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetL(1)+1:end-offsetR(1), offsetL(2)+1:end-offsetR(2), :) = x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :);
% v1 = zeros(sz, c);
v1 = zeros(sz, size(x_, 3));
for l = 1:size(x_, 3)
    v1(:,l) = sdwt2_sara_faceting(x_(:, :, l), I, offset, status, nlevel, wavelet, Ncoefs);
end  
d1 = sqrt(sum(v1.^2,2));
% upsilon = sig*reweight_alpha;
% weights1 = upsilon ./ (upsilon + d1);
weights1 = reweight_alpha ./ (reweight_alpha + d1);
% weights1 = ones(sz, 1);

end
    