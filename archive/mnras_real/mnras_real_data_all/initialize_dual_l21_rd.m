function [v1, weights1] = initialize_dual_l21_rd(x_overlap, ...
    Iq, offset, status_q, nlevel, wavelet, Ncoefs_q, offsetLq, offsetRq, ...
    dims_overlap_ref_q)
% Initialize l21 dual variable for a given facet (constant overlap). Used 
% only for real data (x not zero).
%
% Args:
%     nlevel (int): depth of decomposition.
%     x_overlap (array): overlapping image facet [M, N, L].               
%     Iq (array): starting index of the non-overlapping base facet [1, 2].
%     offset (cell): offset to be used from one dictionary to another 
%                    (different overlap needed for each dictionary -> 
%                    cropping) {nDictionaries}.
%     status_q (array): status of the current facet (last or first 
%                        facet along vert. or hrz. direction) [ndict, 2].
%     nlevel (int): depth of the wavelet decompositions.
%     wavelet (cell): name of the wavelet dictionaries {ndict}.
%     Ncoefs_q (array): size of the wavelet decompositions at each
%                        scale.
%     offsetLq (array): amount of zero-pading from the "left" [1, 2].
%     offsetRq (array): amount of zero-pading from the "right" [1, 2].
%     dims_overlap_ref_q (array): dimension of the facet [1, 2].
%
%
% Returns:
%     v1 (array): dual variable associated with the l21-norm.
%     weights1 (array): weigths ssociated with the l21-norm.
%
%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revised: [../../20..]
%-------------------------------------------------------------------------%
%%

% compute size of the dual variables and the weigths
p = prod(Ncoefs_q, 2);
% if dirac_present
%     sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) + prod(dims);
% else
%     sz = 3*sum(p) - 2*sum(p(nlevel+1:nlevel+1:end));
% end
sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) - 2*p(end); % number of coeffs with the Dirac basis

zerosNum = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap;
v1 = zeros(sz, size(x_, 3));
for l = 1:size(x_, 3)
    v1(:,l) = sdwt2_sara_faceting(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
end  
weights1 = ones(sz, 1);

end
    