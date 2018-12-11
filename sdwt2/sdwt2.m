function [SPsitLx, Ij, dims_PsitLx, Ncoefs] = sdwt2(x_overlap, I, dims, status, lo, hi, J, M)
%
%-------------------------------------------------------------------------%
%%
% Input:
% > x_overlap   overlapping facet extracted from the full image
% > I           start indices of the non-overlapping facets
% > dims        dimension of the underlying non-overlapping facet
% > status      status of the facet considered along each dimension (first: -1, last: 1, none: 0, or both: NaN) )
% > lo          low-pass filter (row vector)
% > hi          high-pass filer (row vector)
% > J           depth of the decomposition
% > M           number of wavelet dictionaries considered (for normalization
%               of the coefficients)
%
% Ouput:
% < SPsitLx     wavelet coefficients [..., 1]
% < Ij          ...
% < dims_PsitLx ...
% < Ncoefs      number of wavelet coefficients associated with the
%               facet of interest (for each level...)
%-------------------------------------------------------------------------%
%%
% Auxiliary variables (sizes)
dim = length(I);
dims_PsitLx = zeros(J+1,2);
Ncoefs = zeros(J+1,2);
Ij = zeros(J+1,2);

% calculating numbers of coefficients in each subband in each dimension [P.-A.] (i.e., h, d, v, a)
m = length(lo);
Sn = I;
Snplus1 = dims + Sn;

for j = 1:J
    for d = 1:dim
        Snj = floor(Sn(d)./2^j);
        if status(d) > 0 || isnan(status(d)) % last / first & last
            Snplus1j = floor(2^(-j).*Snplus1(d)+(1-2^(-j))*(m-1));
        else
            Snplus1j = floor(Snplus1(d)./2^j);
        end
        Ncoefs(j,d) = Snplus1j - Snj; % [P.-A.] (4.19)
        Ij(j,d) = Snj;
    end
end
Ncoefs(J+1,:) = Ncoefs(J,:);
Ij(J+1,:) = Ij(J,:);

for j=1:J-1
    for d=1:dim
        if status(d) < 0 || isnan(status(d)) % first / first & last
            dims_PsitLx(j,d) = Ncoefs(j,d);
        else
            tempDisc = (2^(J-j)-1)*(m-2) + floor(mod(Sn(d),2^J)/2^j); % number of elements to be discarded
            dims_PsitLx(j,d) = Ncoefs(j,d) + tempDisc; % [P.-A.] Nextj (4.21)
        end
    end
end
dims_PsitLx(J+1,:) = Ncoefs(J+1,:);
dims_PsitLx(J,:) = Ncoefs(J,:);
%%%%%%

for i=1:dim
    if dims(i)<2^J
        keyboard
        error('Segment size in dim %d must be >=2^J=%d.',i,2^J);
    end
end

%% forward wavelet transform
PsitLx = cell(3); % (h, v, d) for each 1<= j <= J, a for J+1.
in = x_overlap;
s = 3*sum(prod(Ncoefs(1:end-1,:), 2)) + prod(Ncoefs(end,:)); % total number of coefficients
SPsitLx = zeros(s, 1); % (h, v, d), last row contains only the approximation a

start = 1;
for j=1:J
    % convolution along the rows
    tempa = conv2(in, lo); % conv along rows
    tempd = conv2(in, hi); % extend the signal in a different manner to have another boundary condition: combine wextend and conv(., ., 'valid'), check dimension before that    
    
    % downsampling
    if isnan(status(2)) % first and last
        tempa = tempa(:, 2:2:end);
        tempd = tempd(:, 2:2:end);
    elseif status(2) > 0 % last
        tempa = tempa(:, length(lo):end);
        tempa = tempa(:, 1:2:end);       
        tempd = tempd(:, length(hi):end);
        tempd = tempd(:, 1:2:end); 
    elseif status(2) < 0 % first
        tempa = tempa(:, 1:end-(length(lo)-1));
        tempa = tempa(:, 2:2:end); 
        tempd = tempd(:, 1:end-(length(hi)-1));
        tempd = tempd(:, 2:2:end); 
    else
        tempa = tempa(:, length(lo):end-(length(lo)-1));
        tempa = tempa(:, 1:2:end); 
        tempd = tempd(:, length(hi):end-(length(hi)-1));
        tempd = tempd(:, 1:2:end); 
    end 

    % convolutions along the columns
    PsitLx{1} = conv2(tempa,hi.'); % LH
    PsitLx{2} = conv2(tempd,lo.'); % HL
    PsitLx{3} = conv2(tempd,hi.'); % HH
    in = conv2(tempa,lo.'); % LL
    
    % downsampling
    if isnan(status(1)) % first and last
        PsitLx{1} = PsitLx{1}(2:2:end, :)/sqrt(M);
        PsitLx{2} = PsitLx{2}(2:2:end, :)/sqrt(M);
        PsitLx{3} = PsitLx{3}(2:2:end, :)/sqrt(M);
        in = in(2:2:end, :);
    elseif status(1) > 0 % last
        PsitLx{1} = PsitLx{1}(length(hi):end, :);
        PsitLx{1} = PsitLx{1}(1:2:end, :)/sqrt(M);
        PsitLx{2} = PsitLx{2}(length(lo):end, :);
        PsitLx{2} = PsitLx{2}(1:2:end, :)/sqrt(M);    
        PsitLx{3} = PsitLx{3}(length(hi):end, :);
        PsitLx{3} = PsitLx{3}(1:2:end, :)/sqrt(M);    
        in = in(length(lo):end, :);
        in = in(1:2:end, :);  
    elseif status(1) < 0 % first
        PsitLx{1} = PsitLx{1}(1:end-(length(hi)-1), :);
        PsitLx{1} = PsitLx{1}(2:2:end, :)/sqrt(M); 
        PsitLx{2} = PsitLx{2}(1:end-(length(lo)-1), :);
        PsitLx{2} = PsitLx{2}(2:2:end, :)/sqrt(M); 
        PsitLx{3} = PsitLx{3}(1:end-(length(hi)-1), :);
        PsitLx{3} = PsitLx{3}(2:2:end, :)/sqrt(M); 
        in = in(1:end-(length(lo)-1), :);
        in = in(2:2:end, :); 
    else
        PsitLx{1} = PsitLx{1}(length(hi):end-(length(hi)-1), :);
        PsitLx{1} = PsitLx{1}(1:2:end, :)/sqrt(M); 
        PsitLx{2} = PsitLx{2}(length(lo):end-(length(lo)-1), :);
        PsitLx{2} = PsitLx{2}(1:2:end, :)/sqrt(M); 
        PsitLx{3} = PsitLx{3}(length(hi):end-(length(hi)-1), :);
        PsitLx{3} = PsitLx{3}(1:2:end, :)/sqrt(M); 
        in = in(length(lo):end-(length(lo)-1), :);
        in = in(1:2:end, :); 
    end
    
    % cropping
    sj = prod(Ncoefs(j,:));
    if j < J
        for i = 1:3
            SPsitLx((i-1)*sj+start:i*sj+start-1) = reshape(PsitLx{i}(end-Ncoefs(j,1)+1:end,...
                end-Ncoefs(j,2)+1:end), [sj, 1]);
        end
    else
        for i = 1:3
            SPsitLx((i-1)*sj+start:i*sj+start-1) = reshape(PsitLx{i}, [sj, 1]); % error here...
        end
    end
    start = start + 3*sj;
end

sj = prod(Ncoefs(J+1,:));
SPsitLx(start:start+sj-1) = reshape(in, [sj, 1])/sqrt(M); % approximation

end
