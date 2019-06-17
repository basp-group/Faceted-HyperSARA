function main_scalability_parallel_sdwt2(N, n_runs, n_threads, QyQx)

clc;
format compact
addpath sdwt2

maxNumCompThreads(n_threads);

%% Data and faceting

% Image
% file_name = 'CYGCBEST'; % CYGCBEST
% x = fitsread([file_name,'.fits']);
% x = imresize(x,[512,1024]);   % Resize for Cygnus A
% % x = im2double(x);             % just to make sure that the image is double
% x(x < 0) = 0;
% x(x < 10^(-3)) = 0;
% x = x/max(x(:));  % scale image to avoid numerical errors
% N = size(x);
x = randn(N);

% Algorithm parameters
mu = 1;
gamma_l1 = 1e-3;
rho = 1;
% soft thresholding operator
soft = @(z, T) sign(z) .* max(abs(z)-T, 0);

% Wavelet parameters
n = (1:8);
nlevel = 4;
M = numel(n)+1;
dwtmode('zpd','nodisp');
wavelet = cell(M, 1);
for m = 1:M-1
    wavelet{m} = ['db', num2str(n(m))];
end
wavelet{end} = 'self';
L = [2*n, 0].'; % filter length

%% Time function for a single facet, see time taken (ideally, should be roughly t0/Q)
% Construct Psi for the full image transform
[Psi, Psit] = op_sp_wlt_basis(wavelet, nlevel, N(1), N(2));
Psiu = zeros(N);
u_full = zeros(size(Psit(Psiu)));

% Reference (full image transform)
start_full = tic;
for r = 1:n_runs
    u_old = u_full;
    w = u_full + mu * Psit(x);
    u_full = w - mu * soft(w/mu, gamma_l1/mu) ;
    u_full = rho*u_full + (1-rho)*u_old;
    Psiu = Psi(u_full);
end
t_full = toc(start_full);
t_full = t_full/n_runs;

%% Distributed version

% Create worker pool
% util_create_pool_bis(Q, []);
[~ ,id] = max(prod(QyQx, 2));
Q = prod(QyQx(id, :));
cirrus_cluster = parcluster('local');
cirrus_cluster.NumWorkers = Q;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
if cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads > ncores
    exit(1);
end
% saveProfile(cirrus_cluster);
% clear cirrus_cluster;as
parpool(cirrus_cluster, Q); % override default preference
spmd
    dwtmode('zpd')
end

t_facet = zeros(size(QyQx, 1), 1);

fprintf('================================================================================\n');
fprintf('Faceting primal-dual algorithm \n');
fprintf('================================================================================\n');
fprintf('%10s\t%10s\t%10s\n', 'Qx', 'Qy', 'Acceleration');
fprintf('--------------------------------------------------------------------------------\n');
for k = 1:size(QyQx, 1)
    % Faceting
    Qy = QyQx(k, 1);
    Qx = QyQx(k, 2);
    Q = Qx*Qy;
    rg_y = domain_decomposition(Qy, N(1));
    rg_x = domain_decomposition(Qx, N(2));
    
    segDims = zeros(Q, 4);
    for qx = 1:Qx
        for qy = 1:Qy
            q = (qx-1)*Qy+qy;
            segDims(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1, rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
        end
    end
    I = segDims(:, 1:2);
    dims = segDims(:, 3:4);
    clear segDims
    
    % compute useful indices / sizes
    % Generate all the parameters to define the overlapping facets for each dictionary at each scale
    [~, ~, I_overlap_ref, dims_overlap_ref, I_overlap, dims_overlap, ...
        I_overlap_nc, dims_overlap_nc, status, offset, ~, ~, offsetL, offsetR] = generate_segdwt_indices(N, I, dims, nlevel, wavelet, L);
    u = cell(Q, 1);
    PsiStu = cell(Q, 1);
    for q = 1:Q
        for m = 1:M % to be simplified (precompute sizes to avoid recomputing these later on)
            if ~strcmp(wavelet{m}, 'self')
                [~, Ncoefs_q] = compute_size(I(q, :), dims(q, :), nlevel, status(q,:), L(m)); % not consistent with the size below... to be corrected
                s = 3*sum(prod(Ncoefs_q(1:end-1,:), 2)) + prod(Ncoefs_q(end,:));
                u{q} = [u{q}; zeros(s, 1)];
            else
                u{q} = [u{q}; zeros(prod(dims(q,:)), 1)];
            end
        end
    end
    LtPsiStu = zeros(N);
    
    % Assess scalability of the update step (pseudo-parallel computation)(does not account for communication costs)
    t1 = zeros(n_runs, 1);
    ids = 1:Q;
    for r = 1:n_runs
        it_start = tic;
        updates_L1 = (rand(Q, 1) <= 1); % take all the steps into account
        n_blk = sum(updates_L1);
        idk = ids(updates_L1);
        q1 = 1;
        for q = 1:Q
            if updates_L1(q)
                zerosNum = dims_overlap_ref(q,:) + offsetL(q,:) + offsetR(q,:);
                x_overlap = zeros([zerosNum(:)',1]);
                
                x_overlap(offsetL(q,1)+1:end-offsetR(q,1),...
                    offsetL(q,2)+1:end-offsetR(q,2))...
                    = x(I_overlap_ref(q, 1)+1:I_overlap_ref(q, 1)+dims_overlap_ref(q, 1), ...
                    I_overlap_ref(q, 2)+1:I_overlap_ref(q, 2)+dims_overlap_ref(q, 2));
                
                f(q1) = parfeval(@facet_update, 2, x_overlap, u{q}, I(q,:), dims(q,:), I_overlap_nc{q}, dims_overlap_nc{q}, offset, status(q,:), nlevel, wavelet, N, mu, gamma_l1, rho);
                q1 = q1 + 1;
            end
        end
        
        for q1 = 1:n_blk
            [idx, u_q, PsiStu_q] = fetchNext(f);
            u{idk(idx)} = u_q;
            PsiStu{idk(idx)} = PsiStu_q;
            LtPsiStu = place2DSegment(LtPsiStu, PsiStu{idk(idx)}, I_overlap{idk(idx)}, dims_overlap{idk(idx)});
        end
        t1(r) = toc(it_start);
    end
    t_facet(k) = sum(t1)/n_runs;
    
    fprintf('%10i\t%10i\t%10e\n', Qx, Qy, t_full/t_facet(k));
    fprintf('================================================================================\n');
end

% Assess scalability (parallel version)
save(['results/scalability_parallel_sdwt2_Nx=', num2str(N(2)), '_Ny=', num2str(N(1))], 't_facet', 't_full', 'QyQx')

delete(gcp('nocreate'))
end
