function [PsiSty, I_overlap, dims_overlap] = isdwt2(SPsitLx, I, I_overlap, dims, dims_overlap, Ncoefs, N, lo, hi, J)
%
%-------------------------------------------------------------------------%
%%
% Input:
% > SPsitLx       wavelet coefficients associate dwit the facet considered
% > I             start index of the non-overlapping facet
% > I_overlap     start index of the overlapping facet
% > dims          dimension of the underlying non-overlapping facet
% > dims_overlap  dimension of the overlapping facet
% > Ncoefs        ...
% > N             size of the full image
% > lo            low-pass filter (row vector)
% > hi            high-pass filer (row vector)
% > J             depth of the decomposition
%
% Ouput:
% < PsiSty        
% < I_overlap
% < dims_overlap
%-------------------------------------------------------------------------%
%%
% inverse wavelet transform
dim = 2;
noSubb = dim^2-1;

% precompute rJ
% rJ = zeros(J,dim);
% for j=J-1:-1:1
%     rJ(j,:)=(2^(J-j)-1)*(length(lo)-2) + floor(mod(Sn,2^J)/2^j);
% end
rJ = [(2.^(J-(1:J-1).')-1)*(length(lo)-2) + floor(mod(I,2^J)./(2.^(1:J-1).')); zeros(1, 2)];

% J+1
s = 3*sum(prod(Ncoefs(1:end-1,:), 2)) + prod(Ncoefs(end,:)); % total number of coefficients
sj = prod(Ncoefs(J+1,:));
start = s-sj;
in = reshape(SPsitLx(start+1:sj+start), Ncoefs(J+1,:));

coefTemp = cell(3);
for j=J:-1:2  
    rows = Ncoefs(j-1,1)+rJ(j-1,1);
    sj = prod(Ncoefs(j,:));
    start = start - 3*sj;
    if j == J
        for i=1:noSubb
            coefTemp{i} = reshape(SPsitLx((i-1)*sj+start+1:i*sj+start), Ncoefs(J,:));
        end
        cols = Ncoefs(j,2)+rJ(j,2);
    else
        cols = Ncoefs(j,2)+rJ(j,2);
        for i=1:noSubb
            coefTemp{i} = zeros(rJ(j,:) + Ncoefs(j,:));
            coefTemp{i}(end-Ncoefs(j,1)+1:end, end-Ncoefs(j,2)+1:end) = reshape(SPsitLx((i-1)*sj+start+1:i*sj+start), Ncoefs(j,:));
        end
    end
    
    % upsamling, convolution and cropping along the columns
    tempRows = upConv2c(in, lo, rows, cols) + upConv2c(coefTemp{1}, hi, rows, cols);
    tempRows2 = upConv2c(coefTemp{2}, lo, rows, cols) + upConv2c(coefTemp{3}, hi, rows, cols);
    
    % upsamling, convolution and cropping along the rows
    cols = Ncoefs(j-1,2)+rJ(j-1,2);
    in = upConv2r(tempRows, lo, rows, cols) + upConv2r(tempRows2, hi, rows, cols);  
end


sj = prod(Ncoefs(1,:));
start = start - 3*sj;
for i=1:noSubb
    coefTemp{i} = zeros(rJ(1,:) + Ncoefs(1,:));
    coefTemp{i}(end-Ncoefs(1,1)+1:end, end-Ncoefs(1,2)+1:end) = reshape(SPsitLx((i-1)*sj+start+1:i*sj+start), Ncoefs(1,:));
end

% upsamling, convolution and cropping along the columns
cols = Ncoefs(1,2)+rJ(1,2);
rJ = (2^J-1)*(length(lo)-2)+mod(I(1),2^J);
rows = dims(1)+rJ;
tempRows = upConv2c(in, lo, rows, cols) + upConv2c(coefTemp{1}, hi, rows, cols);
tempRows2 = upConv2c(coefTemp{2}, lo, rows, cols) + upConv2c(coefTemp{3}, hi, rows, cols);

% upsamling, convolution and cropping along the rows
rJ = (2^J-1)*(length(lo)-2) + mod(I(2),2^J);
cols = dims(2) + rJ;
PsiSty = upConv2r(tempRows, lo, rows, cols) + upConv2r(tempRows2, hi, rows, cols); 

%
corner = I_overlap;
dim = length(I_overlap);
dimensions = dims_overlap;

temLIdxs = zeros(dim,1);
temRIdxs = zeros(dim,1);
for i=1:dim
    if corner(i)<0
        temLIdxs(i)= -corner(i);
        dimensions(i) = dimensions(i) - temLIdxs(i);
        corner(i)=0;
    end
    if corner(i)+dimensions(i)>=N(i)
        temRIdxs(i) = corner(i)+dimensions(i) - N(i);
        dimensions(i) = N(i)-corner(i);
    end
end

I_overlap = corner;
dims_overlap = dimensions;

PsiSty = PsiSty(temLIdxs(1)+1:end-temRIdxs(1),temLIdxs(2)+1:end-temRIdxs(2));

end

%% internal function
function y = upConv2c(x, w, rows, cols)

z = zeros(2*size(x, 1), cols);
z(1:2:end, :) = x;
y = conv2(z, w.');
y = y(1:rows, :);

end

function y = upConv2r(x, w, rows, cols)

z = zeros(rows, 2*size(x, 2));
z(:, 1:2:end) = x;
y = conv2(z, w);
y = y(:, 1:cols);

end