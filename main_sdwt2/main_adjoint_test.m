clc; clear all; close all;
format compact

addpath sdwt2

wname = 'db3';

dwtmode('zpd');

N = [256, 125];
nlevel = 5;

X = randn(N);
[PsitX,S] = wavedec2(X,nlevel,wname);
Y = randn(size(PsitX));
PsiY = waverec2(Y,S,wname);

% test adjoint
p1 = sum(conj(PsitX(:).*Y(:)));
p2 = sum(conj(X(:).*PsiY(:)));
abs(p1 - p2) % wavedec and waverec are adjoint for the 'zpd' boundary conditions
             % modifications to waverec are necessary for any other
             % boundary conditions (see codes from Caroline Chaux)
             
%%
% segdwt

% facet definition
% x = double(imread('data/lena.bmp'))/255;
% N = [256, 256];
% N = size(x);
N = [256, 125];
x = randn(N);
Qy = 1;
Qx = 3;
Q = Qx*Qy;
rg_y = domain_decomposition(Qy, N(1));
rg_x = domain_decomposition(Qx, N(2));
randBool = false;

segDims = zeros(Q, 4);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        segDims(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1, rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
    end
end
I = segDims(:, 1:2);
dims = segDims(:, 3:4);

% Wavelet parameters
n = 1:8;
nlevel = 3;
M = numel(n)+1;
dwtmode('zpd','nodisp');
wavelet = cell(M, 1);
for m = 1:M-1
    wavelet{m} = ['db', num2str(n(m))];
end
wavelet{end} = 'self';
L = [2*n,0].'; % filter length

Ij = cell(Q, 1);
dims_PsitLx = cell(Q, 1);
Ncoefs = cell(Q, 1);
% to be possibly simplified: I_overlap_ref_nc/dims_overlap_ref_nc ... Can
% some of these variables be removed
[I_overlap_ref_nc, dims_overlap_ref_nc, I_overlap_ref, dims_overlap_ref, I_overlap, dims_overlap, ...
    I_overlap_nc, dims_overlap_nc, status, offset, ~, ~, offsetL, offsetR] = generate_segdwt_indices(N, I, dims, nlevel, wavelet, L);

u = cell(Q, 1);
for q = 1:Q
    for m = 1:M
        if ~strcmp(wavelet{m}, 'self')
            [~, Ncoefs_q] = compute_size(I(q, :), dims(q, :), nlevel, status(q,:), L(m));
            for j = 1:nlevel
                u{q} = [u{q}; randn(3*prod(Ncoefs_q(j,:)), 1)];
            end
            u{q} = [u{q}; randn(prod(Ncoefs_q(nlevel+1,:)), 1)];
        else
            u{q} = [u{q}; randn(prod(dims(q,:)), 1)];
        end
    end
end

SPsitLx = cell(Q, 1);
PsiStu = cell(Q, 1);
for q = 1:Q
    
    zerosNum = dims_overlap_ref(q,:) + offsetL(q,:) + offsetR(q,:);
    x_overlap = zeros([zerosNum(:)',1]);
    
    x_overlap(offsetL(q,1)+1:end-offsetR(q,1),...
            offsetL(q,2)+1:end-offsetR(q,2))...
            = x(I_overlap_ref(q, 1)+1:I_overlap_ref(q, 1)+dims_overlap_ref(q, 1), ...
        I_overlap_ref(q, 2)+1:I_overlap_ref(q, 2)+dims_overlap_ref(q, 2));
    
    % forward operator [put the following instructions into a parfeval for parallelisation]
    [SPsitLx{q}, Ij{q}, dims_PsitLx{q}, Ncoefs{q}] = sdwt2_sara(x_overlap, I(q, :), dims(q, :), offset, status(q, :), nlevel, wavelet);
    
    % inverse operator (for a single facet) (inverse = adjoin for spd, properly implement the adjoint operator for different boundary conditions) u{q}
    [PsiStu{q}, I_overlap{q}, dims_overlap{q}] = isdwt2_sara(u{q}, I(q, :), dims(q, :), I_overlap_nc{q}, dims_overlap_nc{q}, Ncoefs{q}, N, nlevel, wavelet);
end
% 
LtPsiStu = zeros(N);
for q = 1:Q
    LtPsiStu = place2DSegment(LtPsiStu, PsiStu{q}, I_overlap{q}, dims_overlap{q});
%     imshow(LtPsiStu);
%     pause
end

% norm(LtPsiStu(:) - x(:))

% check adjoint
p1 = sum(conj(x(:).*LtPsiStu(:)));
p2 = 0;

for q = 1:Q
    p2 = p2 + sum(conj(SPsitLx{q}(:)).*u{q}(:));
end

abs(p1 - p2) % ok !

% estimate norm of the operator
val = op_norm_sdwt2(N, I, dims, nlevel, wavelet, L, 1e-4, 100, 0); % val = 1 (as expected)