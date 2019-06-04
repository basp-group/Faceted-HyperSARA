clc; clear all; close all;
format compact;

addpath ../src
addpath ../data
addpath ../utils

% x = double(imread('data/lena.bmp'))/255;
x = fitsread('W28_1024.fits');
N = size(x);

%% Note
% posible errors: value of overlap, value for crop_l21, crop_nuclear, ...
% 

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

%%
% Extract portion of the code to be tested
%%- begin initialization sdwt2
% instantiate auxiliary variables for sdwt2
[~, ~, ~, dims_overlap_ref, I_overlap, dims_overlap, ...
    ~, ~, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices([M, N], I, dims, nlevel, wavelet, L);

% total number of workers (Q: facets workers, K: data workers)
delete(gcp('nocreate'))
numworkers = Q;
cirrus_cluster = parcluster('local');
cirrus_cluster.NumWorkers = numworkers;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
if cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads > ncores
    exit(1);
end

parpool(cirrus_cluster, numworkers);

% define parallel constants (known by each worker)
Qyp = parallel.pool.Constant(Qy);
Qxp = parallel.pool.Constant(Qx);
Qp = parallel.pool.Constant(Q);
waveletp = parallel.pool.Constant(wavelet);
nlevelp = parallel.pool.Constant(nlevel);
offsetp = parallel.pool.Constant(offset);

% define composite variables (local to a given worker)
% wavelet auxiliary variables
Iq = Composite();
dims_q = Composite();
dims_oq = Composite();
temLIdxs_q = Composite();
temRIdxs_q = Composite();
I_overlap_q = Composite();
dims_overlap_q = Composite();
dims_overlap_ref_q = Composite();
status_q = Composite();
Ncoefs_q = Composite();
offsetLq = Composite();
offsetRq = Composite();
% dimension of the ghost cells
overlap_g_south = Composite();
overlap_g_east = Composite();
overlap_g_south_east = Composite();
overlap = Composite();
% 
crop_nuclear = Composite();
crop_l21 = Composite();
%
xq = Comopsite();

% initialize composite variables and constants
for q = 1:Q
    Iq{q} = I(q,:);
    dims_q{q} = dims(q,:);
    dims_oq{q} = dims_o(q,:);
    temLIdxs_q{q} = temLIdxs{q};
    temRIdxs_q{q} = temRIdxs{q};
    I_overlap_q{q} = I_overlap{q};
    dims_overlap_q{q} = dims_overlap{q};
    status_q{q} = status(q,:);
    Ncoefs_q{q} = Ncoefs{q};
    
    % additional composite variables (for the zero padding, see if fewer elements can be used)
    dims_overlap_ref_q{q} = dims_overlap_ref(q,:);
    offsetLq{q} = offsetL(q,:);
    offsetRq{q} = offsetR(q,:);
    
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
    crop_l21{q} = tmp_crop_l21;
    crop_nuclear{q} = tmp_crop_nuclear;    
    overlap{q} = max(max(dims_overlap{q}) - dims(q,:), dims_o(q, :) - dims(q,:));
    
    % initialize the image (portion on each node)
    xq{q} = x(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2));
end

% overlap dimension of the neighbour (necessary to define the ghost cells properly)
w = Composite();
for q = 1:Q
    [qy, qx] = ind2sub([Qy, Qx], q);
    if qy < Qy
        % S (qy+1, qx)
        overlap_g_south{q} = overlap{(qx-1)*Qy + qy+1};
        
        if qx < Qx
            % SE (qy+1, qx+1)
            overlap_g_south_east{q} = overlap{qx*Qy + qy+1};
        else
            overlap_g_south_east{q} = [0, 0];
        end
    else
        overlap_g_south{q} = [0, 0];
        overlap_g_south_east{q} = [0, 0];
    end
    if qx < Qx
        % E (qy, qx+1)
        overlap_g_east{q} = overlap{qx*Qy + qy};
    else
        overlap_g_east{q} = [0, 0];
    end
end

spmd   
    x_overlap = zeros([max_dims, size(xsol_q, 3)]); % zeros([dims_overlap_ref_q, size(xsol_q, 3)]);
    x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = xhat_q;
    x_overlap = comm2d_update_ghost_cells(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
    
    zerosNum = dims_overlap_ref_q + offsetLq + offsetRq;
    x2 = zeros(zerosNum);
    x2(offsetLq(1)+1:end-offsetRq(1),...
            offsetLq(q)+1:end-offsetRq(2), :)...
            = x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :);
    
    % forward operator [put the following instructions into a parfeval for parallelisation]
    SPsitLx = sdwt2_sara(x2, Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q);
    
    % inverse operator (for a single facet) (inverse = adjoin for spd, properly implement the adjoint operator for different boundary conditions) u{q}
    g = isdwt2_sara(SPsitLx, I_q, dims_q, I_overlap_q, dims_overlap_q, Ncoefs_q, nlevelp.Value, waveletp.Value, temLIdxs_q, temRIdxs_q);
    
    % communicate borders + aggregate (reduction)
    g = comm2d_reduce(g, overlap, Qyp.Value, Qxp.Value);
end

x_rec = zeros(N);
for q = 1:Q
    x_rec(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2)) = g{q};
end
figure; imagesc(log10(x_rec));

norm(x(:) - x_rec(:))

