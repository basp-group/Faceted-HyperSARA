function [dims_PsitLx_crop, Ncoefs, Ij] = compute_size(I, dims, J, status, L)
% Compute the number of wavelet coefficients generated for a facet.
%
% Compute the number of wavelet coefficients generated for one facet by the
% faceted wavelet transform described in :cite:`Prusa2012`.
%
% Args:
%     I (array_like): facet start index [1, 2].
%     dims (array_like): dimension of the current facet (w/o overlap) .
%     J (int): number of decomposition levels.
%     status (array_like): status of the facet (first/last) [1, 2] 
%                          value {-1, 0, 1} for (first, none, last).
%     L (int): filter length (wavelet decomposition).
%
% Returns:
%     dims_PsitLx_crop (array_like): dimension of the wavelet coefficients 
%                                    (J+1, 2).
%     Ncoefs (array_like): number of valid coefficients (J+1, 2).
%     Ij (array_like): ... (J+1, 2).
%

%-------------------------------------------------------------------------%
%%
% dim = length(I);
dims_PsitLx_crop = zeros(J+1,2);
Ncoefs = zeros(J+1,2); % dim = 2
Ij = zeros(J+1,2);

% calculating numbers of coefficients in each subband in each dimension [P.-A.] (i.e., h, d, v, a)
Snplus1 = dims + I;

for d = 1:2
    for j = 1:J
        Snj = floor(I(d)./2^j);
        if status(d) > 0 || isnan(status(d)) % last, 1 / first && last
            Snplus1j = floor(2^(-j).*Snplus1(d)+(1-2^(-j))*(L-1));
        else
            Snplus1j = floor(Snplus1(d)./2^j);
        end
        Ncoefs(j,d) = Snplus1j - Snj; % [P.-A.] (4.19)
        Ij(j,d) = Snj;
    end
end
Ncoefs(J+1,:) = Ncoefs(J,:);
Ij(J+1,:) = Ij(J,:);

for d=1:2
    for j=1:J-1
        if status(d) < 0 || isnan(status(d)) % first, -1 / first & last
            dims_PsitLx_crop(j,d) = Ncoefs(j,d);
        else
            tempDisc = (2^(J-j)-1)*(L-2) + floor(mod(I(d),2^J)/2^j); % number of elements to be discarded
            dims_PsitLx_crop(j,d) = Ncoefs(j,d) + tempDisc; % [P.-A.] Nextj (4.21)
        end
    end
end
dims_PsitLx_crop(J+1,:) = Ncoefs(J+1,:);
dims_PsitLx_crop(J,:) = Ncoefs(J,:);

end
