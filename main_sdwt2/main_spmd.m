clc;
format compact

addpath sdwt2

Nx = 4096;
n_runs = 1;
n_threads = 1;
QyQx = [3, 3];
% QyQx = [2, 1; 2, 2; 3, 2; 3, 3; 4, 3; 4, 4];

maxNumCompThreads(n_threads);
N = [Nx, Nx];

%% Data and facetingq

% Image
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

% define parallel constants
waveletp = parallel.pool.Constant(wavelet);
Np = parallel.pool.Constant(N);
gamma_l1p = parallel.pool.Constant(gamma_l1);
rhop = parallel.pool.Constant(rho);
mup = parallel.pool.Constant(mu);
nlevelp = parallel.pool.Constant(nlevel);

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
    rg_y = domain_decomposition2(Qy, N(1));
    rg_x = domain_decomposition2(Qx, N(2));
    
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
    Ij = cell(Q, 1);
    dims_PsitLx = cell(Q, 1);
    Ncoefs = cell(Q, 1);
    [~, ~, I_overlap_ref, dims_overlap_ref, I_overlap, dims_overlap, ...
        I_overlap_nc, dims_overlap_nc, status, offset, ~, ~, offsetL, offsetR] = generate_segdwt_indices(N, I, dims, nlevel, wavelet, L);
    
    % define parallel constant
    offsetp = parallel.pool.Constant(offset);
    
%     u = cell(Q, 1);
    u = Composite(Q);
    x_overlap = Composite(Q);
    PsiStu = Composite(Q);
    Iq = Composite(Q);
    dims_q = Composite(Q);
    I_overlap_nc_q = Composite(Q);
    dims_overlap_nc_q = Composite(Q);
    status_q = Composite(Q);
    for q = 1:Q
        u_q = [];
        for m = 1:M % to be simplified (precompute sizes to avoid recomputing these later on)
            if ~strcmp(wavelet{m}, 'self')
                [~, Ncoefs_q] = compute_size(I(q, :), dims(q, :), nlevel, status(q,:), L(m)); % not consistent with the size below... to be corrected
                s = 3*sum(prod(Ncoefs_q(1:end-1,:), 2)) + prod(Ncoefs_q(end,:));
                u_q = [u_q; zeros(s, 1)];
            else
                u_q = [u_q; zeros(prod(dims(q,:)), 1)];
            end
        end
        u{q} = u_q;
        Iq{q} = I(q, :);
        dims_q{q} = dims(q, :);
        I_overlap_nc_q{q} = I_overlap_nc{q};
        dims_overlap_nc_q{q} = dims_overlap_nc{q};
        status_q{q} = status(q, :);
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
        ticBytes(gcp)
        % test with an spmd instruction
        for q = 1:Q
            zerosNum = dims_overlap_ref(q,:) + offsetL(q,:) + offsetR(q,:);
            x_o = zeros([zerosNum(:)',1]);
            
            x_o(offsetL(q,1)+1:end-offsetR(q,1),...
                offsetL(q,2)+1:end-offsetR(q,2))...
                = x(I_overlap_ref(q, 1)+1:I_overlap_ref(q, 1)+dims_overlap_ref(q, 1), ...
                I_overlap_ref(q, 2)+1:I_overlap_ref(q, 2)+dims_overlap_ref(q, 2));
            x_overlap{q} = x_o;
        end
        tocBytes(gcp)
        
        ticBytes(gcp)
        spmd % why are so many bytes sent here?? (computation supposed to be local only...)
            [PsiStu, u] = facet_update_par(x_overlap, u, Iq, dims_q, I_overlap_nc_q, dims_overlap_nc_q, offsetp, status_q, nlevelp, waveletp, Np, mup, gamma_l1p, rhop);
        end
        tocBytes(gcp)
%         t1(r) = toc(it_start);
        
        ticBytes(gcp)
        for q = 1:Q
            LtPsiStu = place2DSegment(LtPsiStu, PsiStu{q}, I_overlap{q}, dims_overlap{q}); % this part is really slow...
        end
        tocBytes(gcp)
        t1(r) = toc(it_start);
    end
    t_facet(k) = sum(t1)/n_runs;
    
    PsiStu_save = cell(Q, 1);
    PsiStu_save(1:Q) = PsiStu(1:Q);
    
    u_save = cell(Q, 1);
    u_save(1:Q) = u(1:Q);
    
    fprintf('%10i\t%10i\t%10e\n', Qx, Qy, t_full/t_facet(k));
    fprintf('================================================================================\n');
end

% Assess scalability (parallel version)
% save(['results/scalability_parallel_sdwt2_Nx=', num2str(N(2)), '_Ny=', num2str(N(1)), '_test_cst'], 't_facet', 't_full', 'QyQx')

delete(gcp('nocreate'))
