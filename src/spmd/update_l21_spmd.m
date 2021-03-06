function [v1, g] = update_l21_spmd(v1, x_overlap, weights, beta1, Iq, dims_q,...
                                  I_overlap_q, dims_overlap_q, offset, ...
                                  status_q, nlevel, wavelet, Ncoefs_q, ...
                                  temLIdxs_q, temRIdxs_q, offsetLq, ...
                                  offsetRq, dims_overlap_ref_q)
% Update the dual variable associated with the facet l21-norm prior.
%
% Update a facet dual variable associated with the l21-norm prior.
%
% Args:
%     v1 (array_like): dual variable associated with the facet l21-norm 
%                      [s, L].
%     x_overlap (array_like): overlapping image facet [M, N, L].
%     weights (array_like): weights to balance the effect of redundant pixels 
%                           due to the overlap between facets [M, N].
%     beta1 (double): ratio between regularization and convergence 
%                     parameter (gamma1 / sigma1).                   
%     Iq (array_like): starting index of the non-overlapping base facet [1, 2].
%     dims_q (array_like): dimensions of the non-overlapping base facet [1, 2].
%     I_overlap_q (array_like): starting index of the facet [1, 2].
%     dims_overlap_q (array_like): dimensions of the facet [1, 2].
%     offset (cell): offset to be used from one dictionary to another 
%                    (different overlap needed for each dictionary -> 
%                    cropping) {nDictionaries}.
%     status_q (array_like): status of the current facet (last or first 
%                        facet along vert. or hrz. direction) [ndict, 2].
%     nlevel (int): depth of the wavelet decompositions.
%     wavelet (cell): name of the wavelet dictionaries {ndict}.
%     Ncoefs_q ([type]): size of the wavelet decompositions at each
%                        scale.
%     temLIdxs_q (array_like): amount of cropping from the "left" [1, 2].
%     temRIdxs_q (array_like): amount of cropping from the "rights" [1, 2].
%     offsetLq (array_like): amount of zero-pading from the "left" [1, 2].
%     offsetRq (array_like): amount of zero-pading from the "right" [1, 2].
%     dims_overlap_ref_q (array_like): dimension of the facet [1, 2].
%
% Returns:
%     v1 (array_like): dual variable associated with the l21-norm prior 
%                      [s, L].
%     g (array_like): auxiliary variable for the update of the primal 
%                     variable [M, N, L].

%-------------------------------------------------------------------------%
%%
% Code: P.-A. Thouvenin.
% Last revision: [08/08/2019]
%-------------------------------------------------------------------------%
%%

zerosNum = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding
x_ = zeros([zerosNum, size(x_overlap, 3)]);
x_(offsetLq(1)+1:end-offsetRq(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap;
g = zeros(size(x_overlap));

for l = 1 : size(x_, 3)
    w = sdwt2_sara(x_(:, :, l), Iq, offset, status_q, nlevel, wavelet, Ncoefs_q);
    
    w = v1(:, l) +  w;
    l2 = sqrt(sum(abs(w).^2,2));
    l2_soft = max(l2 - beta1*weights, 0)./(l2+eps);
    v1(:, l) = w - l2_soft.*w;
    
    g(:, :, l) = isdwt2_sara(v1(:, l), Iq, dims_q, I_overlap_q, dims_overlap_q, Ncoefs_q, nlevel, wavelet, temLIdxs_q, temRIdxs_q);
end

end
