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
    ~, ~, status, offset, offsetL, offsetR, Ncoefs, temLIdxs, temRIdxs] = generate_segdwt_indices(N, I, dims, 4, {'db8'}, 16);

% % overlap dimension of the neighbour (necessary to define the ghost cells properly)
% overlap = zeros(Q, 2);
% overlap_g_south = zeros(Q, 2);
% overlap_g_east = zeros(Q, 2);
% overlap_g_south_east = zeros(Q, 2);
% for q = 1:Q
%     [qy, qx] = ind2sub([Qy, Qx], q);
%     overlap(q, :) = max(dims_overlap{q}) - dims(q,:);
%     if qy < Qy
%         % S (qy+1, qx)
%         overlap_g_south(q, :) = overlap((qx-1)*Qy + qy+1);
%         
%         if qx < Qx
%             % SE (qy+1, qx+1)
%             overlap_g_south_east(q, :) = overlap(qx*Qy + qy+1);
%         else
%             overlap_g_south_east(q, :) = [0, 0];
%         end
%     else
%         overlap_g_south(q, :) = [0, 0];
%         overlap_g_south_east(q, :) = [0, 0];
%     end
%     if qx < Qx
%         % E (qy, qx+1)
%         overlap_g_east(q, :) = overlap(qx*Qy + qy);
%     else
%         overlap_g_east(q, :) = [0, 0];
%     end
% end

% create weight matrix W (full image size)
W = zeros(N);
for q = 1:Q
   W(I_overlap_ref(q,1)+1:I_overlap_ref(q,1)+dims_overlap_ref(q,1), I_overlap_ref(q,2)+1:I_overlap_ref(q,2)+dims_overlap_ref(q,2)) = ...
       W(I_overlap_ref(q,1)+1:I_overlap_ref(q,1)+dims_overlap_ref(q,1), I_overlap_ref(q,2)+1:I_overlap_ref(q,2)+dims_overlap_ref(q,2)) + ones(dims_overlap_ref(q,:)); 
end
W = 1./W;

figure; imagesc(W); colorbar;

%% run comparison
full_prior = zeros(n_run, 1);
facet_prior = zeros(n_run, 1);
wfacet_prior = zeros(n_run, 1);
facet_prior_wo = zeros(n_run, 1);

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
        facet_prior_wo(n) = facet_prior_wo(n) + facet_nuclear_norm(X(I_overlap_ref(q,1)+1:I_overlap_ref(q,1)+dims(q,1), I_overlap_ref(q,2)+1:I_overlap_ref(q,2)+dims(q,2), :));
    end
end

%% compare results (in average)
m = mean(full_prior);
d1 = mean(abs(full_prior - facet_prior))/m;
d2 = mean(abs(full_prior - wfacet_prior))/m;
d3 = mean(abs(full_prior - facet_prior_wo))/m;
