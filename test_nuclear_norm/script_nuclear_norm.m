clc; clear all; close all
format compact

addpath ../sdwt2


n_run = 1;
rng(0)

N = [1024, 1024];
L = 100;
Qx = 3;
Qy = 3;
Q = Qx*Qy;
nlevel = 4;
filter_size = 16; % use db8 wavelet filter

%% generate auxiliary parameters for the spatial faceting (sdwt2)

% spectral faceting (w/o overlap)
% create starting index of the spatial facets (I) and associated dimensions
% (dims). Beware: I starts at (0, 0)
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
[~, ~, I_overlap_ref, dims_overlap_ref, I_overlap, dims_overlap, ...
    ~, ~, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices(N, I, dims, nlevel, {'db8'}, filter_size);

% create weight matrix W (full image size)
W = zeros(N);
for q = 1:Q
   W(I_overlap_ref(q,1)+1:I_overlap_ref(q,1)+dims_overlap_ref(q,1), I_overlap_ref(q,2)+1:I_overlap_ref(q,2)+dims_overlap_ref(q,2)) = ...
       W(I_overlap_ref(q,1)+1:I_overlap_ref(q,1)+dims_overlap_ref(q,1), I_overlap_ref(q,2)+1:I_overlap_ref(q,2)+dims_overlap_ref(q,2)) + ones(dims_overlap_ref(q,:)); 
end
W = 1./W;

figure; imagesc(W); colorbar;

%% Create overlap with a fixed overlap, lower than or equal to the maximum 
%% overlap size in sdwt2

% maximum overlap in sdwt2, d can be taken smaller in practice
d = (2^nlevel - 1)*(filter_size-1);
rg_yo = domain_decomposition_overlap2(Qy, N(1), d);
rg_xo = domain_decomposition_overlap2(Qx, N(2), d);
Io = zeros(Q, 2);
dims_o = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        Io(q, :) = [rg_yo(qy, 1)-1, rg_xo(qx, 1)-1];
        dims_o(q, :) = [rg_yo(qy,2)-rg_yo(qy,1)+1, rg_xo(qx,2)-rg_xo(qx,1)+1];
    end
end
% create weight matrix Wo (if needed)
Wo = zeros(N);
for q = 1:Q
   Wo(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2)) = ...
       Wo(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2)) + ones(dims_o(q,:)); 
end
Wo = 1./Wo;

figure; imagesc(Wo); colorbar;

%% run comparison
full_prior = zeros(n_run, 1);
facet_prior = zeros(n_run, 1);
wfacet_prior = zeros(n_run, 1);
facet_prior_wo = zeros(n_run, 1);
facet_prior_o = zeros(n_run, 1);
wfacet_prior_o = zeros(n_run, 1);

for n = 1:n_run

    X = randn([N,L]);
    
    %% full prior
    [~,S,~] = svd(reshape(X, [prod(N), L]),'econ');
    full_prior(n) = sum(abs(diag(S)));

    %% faceted prior (overlap from sdwt2)
    facet_prior(n) = 0;
    for q = 1:Q
        facet_prior(n) = facet_prior(n) + facet_nuclear_norm(X(I_overlap_ref(q,1)+1:I_overlap_ref(q,1)+dims_overlap_ref(q,1), I_overlap_ref(q,2)+1:I_overlap_ref(q,2)+dims_overlap_ref(q,2), :));
    end

    %% weighted faceted prior (overlap from sdwt2) -> need overlap size from all the 6 neighbours (not so easy, since varying overlap sizes)
    wfacet_prior(n) = 0;
    for q = 1:Q
        wfacet_prior(n) = wfacet_prior(n) + facet_nuclear_norm( ...
            W(I_overlap_ref(q,1)+1:I_overlap_ref(q,1)+dims_overlap_ref(q,1), I_overlap_ref(q,2)+1:I_overlap_ref(q,2)+dims_overlap_ref(q,2)).* ...
            X(I_overlap_ref(q,1)+1:I_overlap_ref(q,1)+dims_overlap_ref(q,1), I_overlap_ref(q,2)+1:I_overlap_ref(q,2)+dims_overlap_ref(q,2), :));
    end
    
    %% facet prior (no overlap)
    facet_prior_wo(n) = 0;
    for q = 1:Q
        facet_prior_wo(n) = facet_prior_wo(n) + facet_nuclear_norm(X(I(q,1)+1:I(q,1)+dims(q,1), I(q,2)+1:I(q,2)+dims(q,2), :));
    end
    
    %% facet prior (fixed overlap d)
    facet_prior_o(n) = 0;
    for q = 1:Q
        facet_prior_o(n) = facet_prior_o(n) + facet_nuclear_norm(X(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2), :));
    end
    
    %% weighted faceted prior (fixed overlap)
    wfacet_prior_o(n) = 0;
    for q = 1:Q
        wfacet_prior_o(n) = wfacet_prior_o(n) + facet_nuclear_norm( ...
            Wo(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2)).* ...
            X(Io(q,1)+1:Io(q,1)+dims_o(q,1), Io(q,2)+1:Io(q,2)+dims_o(q,2), :));
    end
end

%% compare results (in average)
m = mean(full_prior);
d1 = mean(abs(full_prior - facet_prior))/m;
d2 = mean(abs(full_prior - wfacet_prior))/m;
d3 = mean(abs(full_prior - facet_prior_wo))/m;
d4 = mean(abs(full_prior - facet_prior_o))/m;
d5 = mean(abs(full_prior - wfacet_prior_o))/m;
