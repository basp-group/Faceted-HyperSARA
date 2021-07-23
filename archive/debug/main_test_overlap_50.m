clc; clear all; close all;
format compact;

addpath ../src
addpath ../data
addpath ../utils

x = double(imread('data/lena.bmp'))/255;
% x = fitsread('W28_1024.fits');
N = size(x);

% write a parallel test just in case (spmd) -> later on

%%
% segdwt

% facet definition
Qy = 2;
Qx = 2;
p = 0.5;

% overlapping facets
% define overlapping facets
rg_yo = domain_decomposition_overlap3(Qy, N(1), p);
rg_xo = domain_decomposition_overlap3(Qx, N(2), p);

% equivalent number of facets for the non-redundant version
Qx = floor(Qx/p)+1-floor(1/p);
Qy = floor(Qy/p)+1-floor(1/p);
Q = Qx*Qy;

% define starting indices/sizes of the overlapping facets 
Io = zeros(Q, 2);
dims_o = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        Io(q, :) = [rg_yo(qy, 1)-1, rg_xo(qx, 1)-1];
        dims_o(q, :) = [rg_yo(qy,2)-rg_yo(qy,1)+1, rg_xo(qx,2)-rg_xo(qx,1)+1];
    end
end

% define "base facets"
rg_y = domain_decomposition(Qy, N(1));
rg_x = domain_decomposition(Qx, N(2));
I = zeros(Q, 2);
dims = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        I(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1];
        dims(q, :) = [rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
    end
end

%%
% Wavelet parameters
n = 8;
nlevel = 3;
M = numel(n)+1;
dwtmode('zpd','nodisp');
wavelet = cell(M, 1);
for m = 1:M-1
    wavelet{m} = ['db', num2str(n(m))];
end
wavelet{end} = 'self';
L = [2*n,0].'; % filter length

[I_overlap_ref_nc, dims_overlap_ref_nc, I_overlap_ref, dims_overlap_ref, I_overlap, dims_overlap, ...
    I_overlap_nc, dims_overlap_nc, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices(N, I, dims, nlevel, wavelet, L);

%%
% Extract portion of the code to be tested
crop_l21 = zeros(Q, 2);
crop_nuclear = zeros(Q, 2);
overlap = zeros(Q, 2);
max_dims = zeros(Q, 2);
I_max = zeros(Q, 2);
for q = 1:Q    
    % check which overlap is the largest (sdwt2 or nuclear norms)
    bool_crop = dims_overlap_ref(q,:) >= dims_o(q,:); % 0: nuclear norm largest, 1: dwt2 largest
    tmp_crop_l21 = [0, 0];
    tmp_crop_nuclear = [0, 0];
    if bool_crop(1)
        % sdwt2 largest
        tmp_crop_nuclear(1) = dims_overlap_ref(q,1) - dims_o(q,1);
    else
        tmp_crop_l21(1) = dims_o(q,1) - dims_overlap_ref(q,1);
    end
    if bool_crop(2)
        % sdwt2 largest
        tmp_crop_nuclear(2) = dims_overlap_ref(q,2) - dims_o(q,2);
    else
        tmp_crop_l21(2) = dims_o(q,2) - dims_overlap_ref(q,2);
    end
    crop_l21(q,:) = tmp_crop_l21;
    crop_nuclear(q,:) = tmp_crop_nuclear;    
    overlap(q,:) = max(max(dims_overlap{q}) - dims(q,:), dims_o(q, :) - dims(q,:));   
    max_dims(q,:) = max(dims_overlap_ref(q,:), dims_o(q,:));
    I_max(q,:) = min(I_overlap_ref(q,:), Io(q,:));
end


SPsitLx = cell(Q, 1);
PsiStu = cell(Q, 1);
for q = 1:Q    
%     x_overlap = zeros([max_dims(q,:), size(x, 3)]); % zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
%     x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xq;
    x_overlap = x(I_max(q,1)+1:I_max(q,1)+max_dims(q,1), I_max(q,2)+1:I_max(q,2)+max_dims(q,2), :); % possibly wrong...
    
    zerosNum = dims_overlap_ref(q,:) + offsetL(q,:) + offsetR(q,:); % offset for the zero-padding (to be checked again...)
    x2 = zeros(zerosNum);
    x2(offsetL(q,1)+1:end-offsetR(q,1),...
            offsetL(q,2)+1:end-offsetR(q,2), :)...
            = x_overlap(crop_l21(q,1)+1:end, crop_l21(q,2)+1:end, :);
    
    % forward operator [put the following instructions into a parfeval for parallelisation]
    SPsitLx{q} = sdwt2_sara(x2, I(q, :), offset, status(q, :), nlevel, wavelet, Ncoefs{q}); % check the parallelisation mentioned by Arwa to possibly accelerate the code 
    
    % inverse operator (for a single facet) (inverse = adjoin for spd, properly implement the adjoint operator for different boundary conditions) u{q}
    PsiStu{q} = isdwt2_sara(SPsitLx{q}, I(q, :), dims(q, :), I_overlap{q}, dims_overlap{q}, Ncoefs{q}, nlevel, wavelet, temLIdxs{q}, temRIdxs{q});
end
% 
LtPsiStu = zeros(N);
for q = 1:Q
    %LtPsiStu = place2DSegment(LtPsiStu, PsiStu{q}, min(I_overlap{q}), max(dims_overlap{q})); % check I_overlap{q}, dims_overlap{q}
    LtPsiStu = place2DSegment(LtPsiStu, PsiStu{q}, I_overlap_ref(q, :), dims_overlap_ref(q, :));
    imshow(LtPsiStu);
    pause(1)
end