function PsiSty = isdwt2_sara_new(SPsitLx, I, dims, I_overlap, dims_overlap, Ncoefs, J, wavelet, left_offset, right_offset)
% Inverse operator to compute the contribution of the SARA prior to a 
% single facet.
%-------------------------------------------------------------------------%
%%
% Input:
%
% > SPsitLx          wavelet coefficients obtained from the facet of
%                    interest
% > I                starting index of the facet (without overlap) [1,2]
% > dims             facet dimension (without overlap) [1,2]
% > dims_overlap     dimension of the extended image facets (overlap) [M,2]
% > Ncoefs           dimension of the wavelet facets for each decomposition
%                    level, for each dictionary [M(J+1),2]
%                    {from level 1 to J}
% > J                number of decomposition levels considered
% > wavelet          name of the wavelet dictionary considered {M,1}
% > left_offset      number of coefficients to be cropped from the left of
%                    the reconstructed facet [M,2]
% > right_offset     number of coefficients to be cropped from the right of
%                    the reconstructed facet [M,2]
%
% Output:
%
% < PsiSty  inverse transform
%
%-------------------------------------------------------------------------%
%%
% Debug:
% [23/10/18] ok.
% [19/11/18] code acceleration.
% [18/03/19] further code acceleration [end]
%-------------------------------------------------------------------------%
%%
% Auxiliary variables
M = numel(wavelet);
start = 1;       % current position in the vector SPsitLx (contribution from several wavelet dictionaries)
start_coefs = 1; % current posisiton in the Ncoefs matrix (get number of ceofficients at each scale for a each basis)

% position of the reconstructed dictionary facet inside the global facet
start_facet = I_overlap-min(I_overlap)+1; % [M, 2]
end_facet = start_facet+dims_overlap-1; % [M, 2]

PsiSty = zeros(max(dims_overlap));

for m = 1:M
    % inverse transform
    if ~strcmp(wavelet{m}, 'self')
        [lo_r, hi_r] = wfilters(wavelet{m}, 'r'); % reconstruction filters
        Ncoefs_m = Ncoefs(start_coefs:start_coefs+(J+1)-1,:);
        s = 3*sum(prod(Ncoefs_m(1:end-1, :), 2)) + prod(Ncoefs_m(end,:));
        PsiSty_m = isdwt2(SPsitLx(start:start+s-1), I, ...
            dims, Ncoefs_m, lo_r, hi_r, J, left_offset(m,:), right_offset(m,:));
        start = start + s;
        start_coefs = start_coefs + (J+1);
    else
        s = prod(Ncoefs(start_coefs,:)); % = prod(dims) for the Dirac basis (no boundary extension)
        PsiSty_m = reshape(SPsitLx(start:start+s-1), Ncoefs(start_coefs,:));      
        start = start + s;
        start_coefs = start_coefs + 1;
    end
    % position of the dictionary facet in the larger redundant image facet
    PsiSty(start_facet(m,1):end_facet(m,1), start_facet(m,2):end_facet(m,2)) = ...
    PsiSty(start_facet(m,1):end_facet(m,1), start_facet(m,2):end_facet(m,2)) + PsiSty_m;
end

% Renormalize reconstructed facet (use of several dictionaries)
PsiSty = PsiSty/sqrt(M);

end
