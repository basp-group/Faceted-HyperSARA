function SPsitLx = sdwt2(imApprox, status, wavelet, J, Ncoefs)
%
%-------------------------------------------------------------------------%
%%
% Input:
% > x_overlap   overlapping facet extracted from the full image
% > I           start indices of the non-overlapping facets
% > dims        dimension of the underlying non-overlapping facet
% > status      status of the facet considered along each dimension (first: -1, last: 1, none: 0, or both: NaN) )
% > J           depth of the decomposition
%
% Ouput:
% < SPsitLx     wavelet coefficients [..., 1]
% < Ij          ...
% < dims_PsitLx ...
% < Ncoefs      number of wavelet coefficients associated with the
%               facet of interest (for each level...)
%-------------------------------------------------------------------------%
%% Auxiliary variables (sizes)
[lo, hi] = wfilters(wavelet, 'd'); % decomposition filters
% Ncoefs = zeros(J+1,2);

% calculating numbers of coefficients in each subband in each dimension [P.-A.] (i.e., h, d, v, a)
LoDim = length(lo);
% Snp1   = dims + I;
sBool = ((status > 0) + isnan(status)); % [first and last] or [last]
% for j = 1:J
%     for d = 1:2
%         Snj = floor(I(d)./2^j);
%         if sBool(d) % last / first & last
%             Snp1j = floor(2^(-j).*Snp1(d)+(1-2^(-j))*(LoDim-1)); 
%         else
%             Snp1j = floor(Snp1(d)./2^j);
%         end
%         Ncoefs(j,d) = Snp1j - Snj; % [P.-A.] (4.19)
%     end
% end
% Ncoefs(J+1,:) = Ncoefs(J,:);

%% forward wavelet transform
% total number of coefficients
sJ3     = 3*prod(Ncoefs(1:J,:), 2); %3 * dimensions per level
SPsitLx = zeros(sum(sJ3)+prod(Ncoefs(J+1,:)), 1); % (h, v, d), last row contains only the approximation a
idS1 = LoDim; idS2 = LoDim;
if isnan(status(2)) || status(2) < 0, idS1 = 2; end
if isnan(status(1)) || status(1) < 0, idS2 = 2; end

s = 1;
for j = 1:J
    
    % convolution along the rows
    tempa = conv2(imApprox, lo); % conv along rows
    tempd = conv2(imApprox, hi); % extend the signal in a different manner to have another boundary condition: combine wextend and conv(., ., 'valid'), check dimension before that
    
    % downsampling
    idE1 = size(tempd,2) - (LoDim-1);
    if sBool(2), idE1 = idE1 + (LoDim-1);   end 
    tempa = tempa(:, idS1:2:idE1);
    tempd = tempd(:, idS1:2:idE1);
    
    idE2 = size(tempa,1);
    if sBool(1), idE2 = idE2 + (LoDim-1);  end 
    
    % convolutions along the columns & downsampling
    imApprox = conv2(tempa,lo.');       % LL
    PsitLx = zeros([size(imApprox), 3]);
    PsitLx(:,:,1) = conv2(tempa,hi.'); % LH
    PsitLx(:,:,2) = conv2(tempd,lo.'); % HL
    PsitLx(:,:,3) = conv2(tempd,hi.'); % HH
    PsitLx = PsitLx(idS2:2:idE2, :,:);
    imApprox = imApprox(idS2:2:idE2, :);
    % cropping SPsitLx
    if j < J
        PsitLx = PsitLx(end-Ncoefs(j,1)+1:end, end-Ncoefs(j,2)+1:end,:);
    end
    
    SPsitLx(s:s+sJ3(j)-1) = PsitLx(:);
    s = s + sJ3(j);    
end

SPsitLx(s:end) = imApprox(:);

end
