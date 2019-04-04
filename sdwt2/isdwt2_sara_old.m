function PsiSty = isdwt2_sara(SPsitLx, I, dims, dims_overlap, Ncoefs, J, wavelet, left_offset, right_offset)
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
% Auxiliary vvariables
M = numel(wavelet);
start = 1;       % current position in the vector SPsitLx (contribution from several wavelet dictionaries)
start_coefs = 1; % current posisiton in the Ncoefs matrix (get number of ceofficients at each scale for a each basis)
start_isdwt = 1; % current position in the vector of wavelet coefficients

% Compute total number of coefficients in the reconstruction
id_dirac = find(ismember(wavelet, 'self'), 1);
dirac_present = ~isempty(id_dirac);
if dirac_present
    s = sum(prod(dims_overlap, 2)) - prod(dims_overlap(id_dirac, :)) + prod(dims);
else
    s = sum(prod(dims_overlap, 2));
end
PsiSty = zeros(s, 1);

for m = 1:M
    % inverse transform
    if ~strcmp(wavelet{m}, 'self')
        [lo_r, hi_r] = wfilters(wavelet{m}, 'r'); % reconstruction filters
        Ncoefs_m = Ncoefs(start_coefs:start_coefs+(J+1)-1,:);
        s = 3*sum(prod(Ncoefs_m(1:end-1, :), 2)) + prod(Ncoefs_m(end,:));
        PsiSty_m = isdwt2(SPsitLx(start:start+s-1), I, ...
            dims, Ncoefs_m, lo_r, hi_r, J, left_offset(m,:), right_offset(m,:));
        PsiSty_m = PsiSty_m(:);
        start = start + s;
        start_coefs = start_coefs + (J+1);
    else
        s = prod(Ncoefs(start_coefs,:)); % = prod(dims) for the Dirac basis (no boundary extension)
        PsiSty_m = reshape(SPsitLx(start:start+s-1), [s, 1]);
        start = start + s;
        start_coefs = start_coefs + 1;
    end
    PsiSty(start_isdwt:start_isdwt+numel(PsiSty_m)-1) = PsiSty_m(:);
    start_isdwt = start_isdwt+numel(PsiSty_m);
end

% Renormalize reconstructed facets (use of several dictionaries)
PsiSty = PsiSty/sqrt(M);

end
