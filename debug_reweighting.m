clc; clear all; close all
format compact

addpath src_mnras/
addpath src_mnras/spmd
addpath src_mnras/spmd/weighted
addpath src_mnras/serial
addpath lib/faceted-wavelet-transform/src
addpath lib/faceted-wavelet-transform/lib

rng(1234)
x = randn(1024, 1024, 3);
N = [size(x, 1), size(x, 2)];
nChannel = size(x, 3);

% Wavelet parameters
n = 1:8;
nlevel = 3;
M = numel(n)+1;
wlt_basis = cell(M, 1);
for m = 1:M-1
    wlt_basis{m} = ['db', num2str(n(m))];
end
wlt_basis{end} = 'self';
filter_length = [2*n,0].'; % filter length

% Algorithm parameters
reweight_alpha = 1;
gam = 1e-3;
sigma0 = 1.0;
sigma1 = 1.0;
sigma2 = 1.0;
beta1 = gam/sigma1;
tau = 0.33;
ext_mode = 'zpd';

%% Compare weights obtained with the faceted version, and the original ones
% i.e., compare value of the resulting weighted prior
% see how to recombine the coefficients in the wavelet domain using the different facets
dwtmode(ext_mode, 'nodisp');
[Psi, Psit] = op_sp_wlt_basis(wlt_basis, nlevel, N(1), N(2));
[~, s] = n_wavelet_coefficients(filter_length(1:end-1), N, ext_mode, nlevel);
s = s+prod(N);

global_v1 = zeros(s, nChannel);
for l = 1:nChannel
    global_v1(:, l) = Psit(x(:,:,l));
end

% compute weights
row_norm = sqrt(sum(global_v1.^2, 2));
global_weights = reweight_alpha./(reweight_alpha + row_norm);
[global_v, global_g] = update_dual_l21_serial(global_v1, Psit, Psi, x, global_weights, beta1, sigma1);
global_g = global_g/sigma1;

% r1 = zeros(size(v1));
% g1 = zeros(size(xhat));
% 
% for l = 1:size(xhat, 3)
%     r1(:, l) = v1(:, l) + Psit(xhat(:,:,l));
% end
% l2 = sqrt(sum(abs(r1).^2,2));
% l2_soft = max(l2 - beta1*weights1, 0)./ (l2+eps);
% v1 = r1 - (l2_soft .* r1);
% 
% for l = 1:size(xhat, 3)
%     g1(:,:,l) = sigma1*Psi(v1(:,l));
% end

%% Faceted version

% facets
Qy = 1;
Qx = 2;
Q = Qx*Qy;
rg_y = split_range(Qy, N(1));
rg_x = split_range(Qx, N(2));

I = zeros(Q, 2);
dims = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        I(q,:) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1];
        dims(q,:) = [rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
    end
end

[I_overlap_ref, dims_overlap_ref, I_overlap, dims_overlap, ...
    status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = ...
    sdwt2_setup(N, I, dims, nlevel, wlt_basis, filter_length);

% parallelisation
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
waveletp = parallel.pool.Constant(wlt_basis);
nlevelp = parallel.pool.Constant(nlevel);
offsetp = parallel.pool.Constant(offset);

Iq = Composite();
dims_q = Composite();
temLIdxs_q = Composite();
temRIdxs_q = Composite();
I_overlap_q = Composite();
dims_overlap_q = Composite();
dims_overlap_ref_q = Composite();
status_q = Composite();
Ncoefs_q = Composite();
offsetLq = Composite();
offsetRq = Composite();
overlap_g_south = Composite();
overlap_g_east = Composite();
overlap_g_south_east = Composite();
overlap = Composite();
crop_l21 = Composite();
for q = 1:Q
    Iq{q} = I(q, :);
    dims_q{q} = dims(q, :);
    temLIdxs_q{q} = temLIdxs{q};
    temRIdxs_q{q} = temRIdxs{q};
    I_overlap_q{q} = I_overlap{q};
    dims_overlap_q{q} = dims_overlap{q};
    status_q{q} = status(q, :);
    Ncoefs_q{q} = Ncoefs{q};
    
    % additional composite variables (for the zero padding / boundary conditions)
    dims_overlap_ref_q{q} = dims_overlap_ref(q,:);
    offsetLq{q} = offsetL(q,:);
    offsetRq{q} = offsetR(q,:);
    crop_l21{q} = [0, 0];  
    overlap{q} = max(dims_overlap{q}) - dims(q,:); 
end

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

%%
x_q = Composite();
g_q = Composite();
for q = 1:Q
    x_q{q} = x(I(q, 1)+1:I(q, 1)+dims(q, 1), I(q, 2)+1:I(q, 2)+dims(q, 2), :);
    g_q{q} = zeros([dims(q, :), nChannel]);
end

% Update constant dual variables
beta1_ = parallel.pool.Constant(gam/sigma1);

id_dirac = find(ismember(wlt_basis, 'self'), 1);
dirac_present = ~isempty(id_dirac);

%%
% Primal / prior nodes (l21/nuclear norm dual variables)
reweight_alphap = Composite();
for q = 1:Q
    reweight_alphap{q} = reweight_alpha;
end
v1_ = Composite();
weights1_ = Composite();
spmd
    dwtmode(ext_mode,'nodisp');
    max_dims = dims_overlap_ref_q;
    p = prod(Ncoefs_q, 2);
%     if dirac_present
%         sz = 3*sum(p(1:end)) -
%         2*sum(p(nlevelp.Value+1:nlevelp.Value+1:end)) + prod(dims_q); %
%         wrong
%     else
%         sz = 3*sum(p) - 2*sum(p(nlevelp.Value+1:nlevel+1:end));
%     end
    sz = 3*sum(p(1:end)) - 2*sum(p(nlevel+1:nlevel+1:end)) - 2*p(end);
    
    % update weights
    x_overlap = zeros([max_dims, size(x_q, 3)]);
    x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = x_q;
    x_overlap = comm2d_update_borders(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
    
    % initialize v1_ and weights1_
    zerosNum = dims_overlap_ref_q + offsetLq + offsetRq; % offset for the zero-padding
    x_ = zeros([zerosNum, size(x_q, 3)]);
    x_(offsetLq(1)+1:end-offsetR(1), offsetLq(2)+1:end-offsetRq(2), :) = x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :);
    v1_ = zeros(sz, size(x_, 3));
    for l = 1 : size(x_, 3)
        v1_(:, l) = sdwt2_sara_faceting(x_(:, :, l), Iq, offsetp.Value, status_q, nlevelp.Value, waveletp.Value, Ncoefs_q);
    end
    d1 = sqrt(sum(v1_.^2,2));
    weights1_ = reweight_alpha ./ (reweight_alpha + d1);

end

%%

spmd
    x_overlap(overlap(1)+1:end, overlap(2)+1:end, :) = x_q;
    x_overlap = comm2d_update_borders(x_overlap, overlap, overlap_g_south_east, overlap_g_south, overlap_g_east, Qyp.Value, Qxp.Value);
    [v1_, g1] = update_dual_l21(v1_, x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :), weights1_, beta1_.Value, Iq, ...
                dims_q, I_overlap_q, dims_overlap_q, offsetp.Value, status_q, ...
                nlevelp.Value, waveletp.Value, Ncoefs_q, temLIdxs_q, temRIdxs_q, offsetLq, offsetRq, dims_overlap_ref_q);
    g = zeros(size(x_overlap));
    g(crop_l21(1)+1:end, crop_l21(2)+1:end, :) = g1(crop_l21(1)+1:end, crop_l21(2)+1:end, :); % x_overlap(crop_l21(1)+1:end, crop_l21(2)+1:end, :);
    g_q = g(overlap(1)+1:end, overlap(2)+1:end, :); % no overlap here 
end

%%
% gather all image g segments in a single image
faceted_g = zeros(size(x));
for q = 1:Q
    faceted_g = place2DSegment(faceted_g, g{q}, min(I_overlap{q},[],1), max(dims_overlap{q},[],1));
end

%% check error between the faceted and serial update step for the l21 variable 
% after a single reweighting step

err = norm(faceted_g(:) - global_g(:));
