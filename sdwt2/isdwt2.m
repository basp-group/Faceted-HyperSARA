function PsiSty = isdwt2(SPsitLx, I, dims, Ncoefs, wavelet, J, temLIdxs, temRIdxs)
%% inverse wavelet transform
[lo, hi] = wfilters(wavelet, 'r');
SzdJ = Ncoefs(1:J, :) + [(2.^(J-(1:J-1).')-1)*(length(lo)-2) + floor(mod(I, 2^J)./(2.^(1:J-1).')); zeros(1, 2)];
dS3  = 3*prod(Ncoefs(1:J, :), 2);

% case J+1
Pos = sum(dS3);
PsiSty = reshape(SPsitLx(1+Pos : prod(Ncoefs(J+1, :))+Pos), Ncoefs(J+1, :));

% case J
dSPsitLx = reshape(SPsitLx(Pos - dS3(J)+1 : Pos),[Ncoefs(J,:), 3]);
Pos  = Pos - dS3(J);
idx = SzdJ(J-1, :);
% upsampling, convolution and cropping
PsiSty = upConv2r(upConv2c(PsiSty, lo, idx(1)) + upConv2c(dSPsitLx(:, :, 1), hi, idx(1)), lo, idx(2)) + ...
    upConv2r(upConv2c(dSPsitLx(:, :, 2), lo, idx(1)) + upConv2c(dSPsitLx(:, :, 3), hi, idx(1)) ,hi, idx(2));

% case 1<j<J
for j = J-1:-1:2
    dSPsitLx = zeros([SzdJ(j, :), 3]);
    dSPsitLx(end-Ncoefs(j,1)+1:end, end-Ncoefs(j,2)+1:end, :) = ...
        reshape(SPsitLx(Pos - dS3(j)+1 : Pos),[Ncoefs(j,:), 3]);
    Pos  = Pos - dS3(j);
    idx = SzdJ(j-1, :);
    % upsampling, convolution and cropping
    PsiSty = upConv2r(upConv2c(PsiSty, lo, idx(1)) + upConv2c(dSPsitLx(:,:,1), hi, idx(1)), lo, idx(2)) + ...
        upConv2r(upConv2c(dSPsitLx(:, :, 2), lo, idx(1)) + upConv2c(dSPsitLx(:,:,3), hi, idx(1)), hi, idx(2));
end

% j = 1
dSPsitLx = zeros([SzdJ(1,:), 3]);
dSPsitLx(end-Ncoefs(1,1)+1:end, end-Ncoefs(1,2)+1:end, :) = ...
    reshape(SPsitLx(Pos - dS3(1)+1 : Pos),[Ncoefs(1,:), 3]);
idx = dims + (2^J-1)*(length(lo)-2) + mod(I,2^J);
% upsampling, convolution and cropping    
PsiSty = upConv2r(upConv2c(PsiSty, lo, idx(1)) + upConv2c(dSPsitLx(:, :, 1), hi, idx(1)), lo, idx(2)) + ...
    upConv2r(upConv2c(dSPsitLx(:, :, 2), lo, idx(1)) + upConv2c(dSPsitLx(:, :, 3), hi, idx(1)), hi, idx(2));  

%% Crop
PsiSty = PsiSty(temLIdxs(1):size(PsiSty,1)-temRIdxs(1),temLIdxs(2):size(PsiSty,2)-temRIdxs(2));
end

%% internal function
function y = upConv2c(x, w, rows)
z = zeros(2*size(x, 1), size(x,2));
z(1:2:end, :) = x;
y = conv2(z, w.');
y = y(1:rows, :);
end

function y = upConv2r(x, w, cols)
z = zeros(size(x,1), 2*size(x, 2));
z(:, 1:2:end) = x;
y = conv2(z, w);
y = y(:, 1:cols);

end