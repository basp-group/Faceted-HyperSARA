function [v0, v1, weights0, weights1] = initialize_dual_and_weights(x_overlap, ...
    I, offset, status, nlevel, wavelet, Ncoefs, dims_o, c, dims_overlap_ref, ...
    offsetL, offsetR, reweight_alpha, crop_l21, crop_nuclear, w, sig, sig_bar)
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

%! dual variables initialized to 0, and weights initialized from the 
%! current value of the primal variable

%! checking if x_overlap = 0 is not strictly speaking (weights are automatically
%! 1 in this case), but it can save some operations

% compute size of the dual variables and the weigths
p = prod(Ncoefs, 2);
% if dirac_present
%     sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) + prod(dims);
% else
%     sz = 3*sum(p) - 2*sum(p(nlevel+1:nlevel+1:end));
% end
% number of wavelet coeffs when the Dirac basis is used
sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) - 2*p(end);
flag_zero = (norm(x_overlap(:)) == 0);

% nuclear norm
v0 = zeros(prod(dims_o), c);
if flag_zero
    weights0 = ones(min(prod(dims_o), c), 1);
else
    sol = w.*x_overlap(crop_nuclear(1)+1:end, crop_nuclear(2)+1:end, :);
    sol = reshape(sol, [numel(sol)/size(sol, 3), size(x_overlap, 3)]);
    [~,S00,~] = svd(sol,'econ');
    d0 = abs(diag(S00));
    upsilon_bar = sig_bar*reweight_alpha;
    weights0 = upsilon_bar ./ (upsilon_bar + d0);
end

% l21 norm
v1 = zeros(sz, c);
if flag_zero
    weights1 = ones(sz, 1);
else
    zerosNum = dims_overlap_ref + offsetL + offsetR; % offset for the zero-padding
    x_ = zeros([zerosNum, size(x_overlap, 3)]);
    x_(offsetL(1)+1:end-offsetR(1), offsetL(2)+1:end-offsetR(2), :) = x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :);
    z = zeros(sz, c);
    for l = 1 : size(x_, 3)
        z(:, l) = sdwt2_sara_faceting(x_(:, :, l), I, offset, status, nlevel, wavelet, Ncoefs);
    end
    d1 = sqrt(sum(z.^2,2));
    upsilon = sig*reweight_alpha;
    weights1 = upsilon ./ (upsilon + d1);
end

end
    