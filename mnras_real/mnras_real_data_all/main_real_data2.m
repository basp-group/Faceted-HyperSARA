function main_real_data2(Qx, Qy, Qc, ...
    algo_version, window_type, ncores_data, ind, overlap_size, ...
    nReweights, extract_raw_data, flag_computeOperatorNorm, ...
    flag_solveMinimization, gam0, gam, rw, alpha, flag_nonZeroPrimal, flag_homotopy)
% Main script to run the faceted HyperSARA approach on synthetic data.
% 
% This script generates synthetic data and runs the faceted HyperSARA 
% approach to reconstruct an :math:`N \times L` wideband image cube.
% 
% Args:
%     Qx (int): number of spatial facets along axis 2 (x)
%     Qy (int): number of spatial facets along axis 1 (y)
%     Qc (int): number of spectral facets
%     algo_version (string): selected version of the solver:
%        - 's'             'standard': overlap of the faceted nuclear norm 
%                          equal to the
%                          one needed for the faceted SARA dictionary, w/o 
%                          spatial weights (apodization window)
%        - 's2'            'standard2': alternative parallellization scheme 
%                          compared to 'standard' (gather image and data on 
%                          the same nodes)
%        - 'c'             'constant': constant overlap for the faceted 
%                          nuclear norm, w/o spatial weights (apodization 
%                          window)
%        - 'w'             'weighted': overlap of the faceted nuclear norm 
%                          equal to the one needed for the faceted SARA 
%                          dictionary,  using spatial weights (apodization 
%                          window)
%        - 'cw'            cst_weighted: constant overlap taken for the 
%                          faceted nuclear norm, using spatial weights 
%                          (apodization window)
%        - 'no'            no_overlap: no overlap for the faceted nuclear 
%                          norm prior
%     window_type (string): type of apodization window considered for the 
%                           faceted nuclear norm prior. Only active with 
%                           the following versions of the algorithm:
%        - 'triangular'
%        - 'hamming'
%        - 'pc' (piecewise-constant)
%     ncores_data (int): number of cores handlig the data fidelity terms 
%                        ("data cores"). The total number of cores is 
%                        Qx*Qy + ncores_data + 1
%     ind (int): index of the spectral facet to be reconstructed (set to -1
%                to deactivate spectral faceting)
%     overlap_size (int): number of overlapping pixels between contiguous 
%                         facets (only active for  the 'cst' and 
%                         'cst_overlap' versions of the solver)
%     nReweights (int): number of reweighting steps
%     flag_extractData (bool): flag specifying whether the real data need 
%       to be extracted
%     flag_computeOperatorNorm (bool): compute norm of the measurement 
%     operator, or load it from an existing .mat file (e.g., computed from 
%     a previous run)
%     flag_solveMinimization (bool): flag triggering the faceted HyperSARA 
%     solver
%     gam (double): l21-norm regulariation parameter
%     gam0 (double): nuclear-norm regularization parameter
%
%%
format compact;

rootgroup = settings();
rootgroup.matlab.general.matfile.SaveFormat.PersonalValue = 'v7.3';

disp('MNRAS (real data) configuration')
disp(['Algorithm version: ', algo_version]);
disp(['Number of facets Qy x Qx : ', num2str(Qy), ' x ', num2str(Qx)]);
disp(['Number of spectral facets Qc : ', num2str(Qc)]);
disp(['Number of cores data fidelity : ', num2str(ncores_data)]);
disp(['Overlap size: ', num2str(overlap_size)]);
disp(['Extracting real data: ', num2str(extract_raw_data)]);
disp(['Computing operator norm: ', num2str(flag_computeOperatorNorm)]);
disp(['Solving problem: ', num2str(flag_solveMinimization)]);

addpath ../../lib/generate_data/
addpath ../../lib/operators/
addpath ../../lib/measurement-operator/nufft/
addpath ../../lib/utils/
addpath ../../lib/faceted-wavelet-transform/src
addpath ../../data/
addpath ../../src_mnras/
% only cw version used in this script
addpath ../../src_mnras/spmd
addpath ../../src_mnras/spmd/weighted

%! TO BE CHANGED (paths, flags, ...)
% setting paths to results and reference image cube
generate_eps_nnls = false;
save_data = false; 
save_full_operator = false;
%extract_raw_data = false;
data_path = '../../data/';
results_path = fullfile('results/');
mkdir(data_path)
mkdir(results_path)
%!----

%% Default config parameters (can be kept as is)
% gridding parameters
Nx = 2560;
Ny = 1536;
N = Nx * Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 7; % number of neighbours for nufft
Ky = 7; % number of neighbours for nufft

% preconditioning parameters
param_precond.N = N;       % number of pixels in the image
param_precond.Nox = ox*Nx; % number of Fourier points (oversampled plane)
param_precond.Noy = oy*Ny;
param_precond.gen_uniform_weight_matrix = 1; % set weighting type
param_precond.uniform_weight_sub_pixels = 1;

% block structure
param_block_structure.use_density_partitioning = 0;
param_block_structure.density_partitioning_no = 1;
param_block_structure.use_uniform_partitioning = 0;
param_block_structure.uniform_partitioning_no = 4;
param_block_structure.use_equal_partitioning = 0;
param_block_structure.equal_partitioning_no = 1;
param_block_structure.use_manual_frequency_partitioning = 0;
% sparam.fpartition = [pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [0, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.25*pi, 0, 0.25*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-64/256*pi, 0, 64/256*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi/4), 0, icdf('norm', 0.75, 0, pi/4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_block_structure.use_manual_partitioning = 1;

%% Generate/load real data
%! TO BE CHECKED / CHANGED
%if flag_extractData
    % visibility_file_name = '/home/h019/aa16/S-band/CYG-S.mat';
    visibility_file_name = {'5','7'};
    configuration_file_name = {'C','A'};

    % load real_data;
    if extract_raw_data
        [y, uw1, vw1, nWw, f, time_, pos] = util_load_real_vis_data(visibility_file_name,configuration_file_name, ind);
        % y
        % norm(y{1})
        % norm(y{32})
        [uw, vw, pixel_size] = scale_uv(visibility_file_name,configuration_file_name, uw1, vw1);
        save(['CYG_data_raw_ind=' num2str(ind) '.mat'],'-v7.3', 'y', 'uw', 'vw', 'nWw', 'f', 'time_', 'pos', 'pixel_size');
    else
        load(['CYG_data_raw_ind=' num2str(ind) '.mat']);
    end

    % figure(1),
    % for i = length(uw)
    %     scatter(uw{i},vw{i},'.');hold on;
    % end
    %ch = [1 : size(y,2)];
    %ch = 31;
    ch = [1:17, 19, 21:32];
    ii = 1;

    for i = ch
        i
        
        %% compute weights
        [aWw] = util_gen_preconditioning_matrix(uw{i}, vw{i}, param_precond);
        
        % set the blocks structure
        param_block_structure.partition = pos{i};
        [u, v, ~, uvidx, aW{ii}, nW] = util_gen_block_structure(uw{i}, vw{i}, aWw, nWw{i}, param_block_structure);
        u
        v
        % measurement operator initialization
        fprintf('Initializing the NUFFT operator\n\n');
        tstart = tic;
        
        % compute A & At with Kx = Ky = 7 
        [A, At, ~, ~] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW); % Kx = Ky = 7
        
        % load the G matrices and the estimated vis.
        if i < 17
            load(['/lustre/home/shared/sc004/AD_Update_G_Codes/results/SubCube' num2str(ind) '/WB5-' num2str(i) '/PreProcStruct.mat']);
        else
            load(['/lustre/home/shared/sc004/AD_Update_G_Codes/results/SubCube' num2str(ind) '/WB7-' num2str(i) '/PreProcStruct.mat']);
        end
        PreProcStruct.Gw
        Gw = [PreProcStruct.Gw{2}; PreProcStruct.Gw{1}]; % concatenate the configurations C first and then A
        res_est_full = [PreProcStruct.ResiudalModelDataInit{2}; PreProcStruct.ResiudalModelDataInit{1}]; % concatenate the configurations C first and then A
        param_HSI.xsol_Arwa(:,:,ii) = PreProcStruct.SolInit;
        clear PreProcStruct;
        disp('Gw')
        size(Gw)

        G{ii} = cell(length(u),1);
        yb{ii} = cell(length(u),1);
        res_est{ii} = cell(length(u),1);
        eps_b{ii} = cell(length(u),1);

        W{ii}  =cell(length(u),1); %%%%%%%% ARWA
        for j = 1 : length(u)
            G{ii}{j} = Gw(uvidx{j},:);
            W{ii}{j}= G{ii}{j}' * ones(size(G{ii}{j}, 1), 1) ~= 0;
            G{ii}{j} = G{ii}{j}(:, W{ii}{j});
            res_est{ii}{j} = res_est_full(uvidx{j});
            yb{ii}{j} = y{i}(uvidx{j});
            eps_b{ii}{j} = norm(res_est{ii}{j});
        end
        clear Gw res_est_full;

        ii = ii + 1;
    end
    save(['eps_new_ind=' num2str(ind) '.mat'],'-v7.3', 'eps_b');

    %% Save data
    if save_data
        save('CYG_data.mat','-v7.3', 'G', 'W', 'aW', 'yb');
        save('CYG_y.mat','-v7.3', 'y');
    end

    %%
    if save_full_operator && exist('Gw','var')
        save('CYG_Gw.mat','-v7.3', 'Gw');
    end

    %% Free memory
    clear y u v uv_mat uv_mat1 uvidx uw vw ant1 ant2 aWw nW nWw out_block;
    clear param_block_structure param_precond;
%else
%    load('CYG_data.mat');
%    [A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);
%end

ch = (1 : size(yb,2));
nChannels = size(yb,2)
%!-----------

%% Setup name of results file
%! TO BE CHANGED
data_name_function = @(nchannels) strcat('y_cygA_Nx=',num2str(Nx),'_Ny=',num2str(Ny), '_L=', num2str(nchannels),'.mat');

results_name_function = @(nchannels) strcat('Final_fhs_', algo_version,'_',window_type,'_Nx=',num2str(Nx), '_Ny=',num2str(Ny), ...
    '_L=',num2str(nchannels), '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), ...
    '_Qc=', num2str(Qc), '_overlap=', num2str(overlap_size), '_ind=', num2str(ind), ...
    '_gam0=', num2str(gam0), '_gam=', num2str(gam), '_rw=', num2str(rw), ...
    '_primal=', num2str(flag_nonZeroPrimal), '_homotopy=', num2str(flag_homotopy), '.mat');

temp_results_name = @(nchannels) strcat('fhs_', algo_version,'_',window_type,'_Nx=',num2str(Nx), '_Ny=',num2str(Ny), ...
    '_L=',num2str(nchannels), '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), ...
    '_Qc=', num2str(Qc), '_overlap=', num2str(overlap_size), '_ind=', num2str(ind), ...
    '_gam0=', num2str(gam0), '_gam=', num2str(gam), ...
    '_primal=', num2str(flag_nonZeroPrimal), '_homotopy=', num2str(flag_homotopy));

warm_start = @(nchannels) strcat('fhs_', algo_version,'_',window_type,'_Nx=',num2str(Nx), '_Ny=',num2str(Ny), ...
    '_L=',num2str(nchannels), '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), ...
    '_Qc=', num2str(Qc), '_overlap=', num2str(overlap_size), '_ind=', num2str(ind), ...
    '_gam0=', num2str(gam0), '_gam=', num2str(gam), '_rw=', num2str(rw), ...
    '_primal=', num2str(flag_nonZeroPrimal), '_homotopy=', num2str(flag_homotopy), '.mat');

data_name = data_name_function(nChannels);    
results_name = results_name_function(nChannels);
%!-----------


%% Compute operator norm
%! TO BE CHANGED (paths, ...) add computation of the operator norm w/o 
%! preconditioning?
if flag_computeOperatorNorm
    % Compute full measurement operator spectral norm
    F = afclean( @(x) HS_forward_operator_precond_G(x, G, W, A, aW));
    Ft = afclean( @(y) HS_adjoint_operator_precond_G(y, G, W, At, aW, Ny, Nx));
    Anorm = pow_method_op(F, Ft, [Ny Nx nchans]);
    %save(fullfile(results_path,strcat('Anorm_N=',num2str(Nx), ...
    %    '_L=',num2str(nChannels),'_Qc=',num2str(Qc),'_ind=',num2str(ind), '.mat')),'-v7.3', 'Anorm');
    % save(['Anorm_ind=' num2str(ind) '.mat'],'-v7.3', 'Anorm');
else
    %load(fullfile(results_path,strcat('Anorm_N=',num2str(Nx),...
    %    '_L=',num2str(nChannels),'_Qc=',num2str(Qc),'_ind=',num2str(ind), '.mat')));
    load(['Anorm_ind=' num2str(ind) '.mat']); 
end
%!----------

%% Generate initial eps_b by performing imaging with NNLS on each data block separately
if generate_eps_nnls
    % param_nnls.im = im; % original image, used to compute the SNR
    param_nnls.verbose = 2; % print log or not
    param_nnls.rel_obj = 1e-5; % stopping criterion
    param_nnls.max_iter = 1000; % max number of iterations
    param_nnls.sol_steps = [inf]; % saves images at the given iterations
    param_nnls.beta = 1;
    % solve nnls per block
    parfor i = ch
        eps_b{i} = cell(length(G{i}),1);
        for j = 1 : length(G{i})
            % printf('solving for band %i\n\n',i)
            [~,eps_b{i}{j}] = fb_nnls_blocks(yb{i}{j}, A, At, G{i}{j}, W{i}{j}, param_nnls);
        end
    end
    %! TO BE UPDATED
    save(['eps_ind=' num2str(ind) '.mat'],'-v7.3', 'eps_b');
    %! TO BE UPDATED
    %load(['eps_ind=' num2str(ind) '.mat']);
end

%% Solver
if flag_solveMinimization
    % Definition of the SARA dictionary
    nlevel = 4; % depth of the wavelet decompositions
    %! always specify Dirac basis ('self') in last position if used in the
    %! SARA dictionary 
    wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; 
    filter_length = [2*(1:8)'; 0]; % length of the filters (0 corresponding to the 'self' basis)

    % compute sig and sig_bar (estimate of the "noise level" in "SVD" and 
    % SARA space) involved in the reweighting scheme
    [sig, sig_bar, max_psf, ~, ~, ~] = compute_reweighting_lower_bound(yb, W, G, A, At, Ny, Nx, oy, ox, ...
    nChannels, wlt_basis, filter_length, nlevel);
    % mu = nuclear_norm/l21_norm;
    
    %% HSI parameter structure sent to the  HSI algorithm
    param_HSI.verbose = 2; % print log or not
    param_HSI.nu0 = 1; % bound on the norm of the Identity operator
    param_HSI.nu1 = 1; % bound on the norm of the operator Psi
    param_HSI.gamma0 = gam0; % regularization parameter nuclear norm
    param_HSI.gamma = gam;  % regularization parameter L21 (soft th parameter)
    param_HSI.rel_var = 1e-6;  % stopping criterion
    param_HSI.max_iter = 100000; % max number of iterations
    
    param_HSI.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
    param_HSI.adapt_eps_start = 500; % minimum num of iter before stating adjustment
    param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
    param_HSI.adapt_eps_tol_out = 1.01; % tolerance outside the l2 ball
    param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
    param_HSI.adapt_eps_rel_var = 5e-4; % bound on the relative change of the solution
    param_HSI.adapt_eps_change_percentage = (sqrt(5)-1)/2; % the weight of the update w.r.t the l2 norm of the residual data
    
    %! TO BE CHECKED: see where reweight_alpha needs to start from
    if flag_homotopy
        param_HSI.reweight_alpha = 10;
        param_HSI.reweight_alpha_ff = (1/param_HSI.reweight_alpha)^(1/6); % reach the floor level in 6 reweights (see if a different number would be appropriate)
        % param_HSI.reweight_alpha_ff = 0.8;
        % param_HSI.reweight_alpha = (0.8)^alpha; % 1 % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
    else
        param_HSI.reweight_alpha = 1;
        param_HSI.reweight_alpha_ff = 1;
    end
    param_HSI.total_reweights = nReweights; % -1 if you don't want reweighting
    param_HSI.reweight_abs_of_max = Inf; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised
    param_HSI.sig = sig; % estimate of the noise level in SARA space
    param_HSI.sig_bar = sig_bar; % estimate of the noise level in "SVD" space
    param_HSI.max_psf = max_psf;

    param_HSI.use_reweight_steps = 0; % reweighting by fixed steps
    param_HSI.reweight_step_size = 300; % reweighting step size
    param_HSI.reweight_steps = (5000: param_HSI.reweight_step_size :10000);
    param_HSI.step_flag = 1;

    param_HSI.use_reweight_eps = 1; % reweighting w.r.t the relative change of the solution
    param_HSI.reweight_max_reweight_itr = param_HSI.max_iter - param_HSI.reweight_step_size;
    param_HSI.reweight_rel_var = 5e-4; % criterion for performing reweighting
    param_HSI.reweight_min_steps_rel_obj = 300; % min num of iter between reweights
    
    param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid   
  
    %% faceted HyperSARA
    param_HSI.nu2 = Anorm; % upper bound on the norm of the measurement operator A*G
    param_HSI.num_workers = Qx*Qy + ncores_data;

    disp('Faceted HyperSARA')
    disp('-----------------------------------------')

    %! spectral tesselation (non-overlapping)
    rg_c = split_range(ncores_data, ch(end));
    cell_c_chunks = cell(ncores_data, 1);
    y_spmd = cell(ncores_data, 1);
    epsilon_spmd = cell(ncores_data, 1);
    aW_spmd = cell(ncores_data, 1);
    W_spmd = cell(ncores_data, 1);
    G_spmd = cell(ncores_data, 1);

    for i = 1:ncores_data
        cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
        y_spmd{i} = yb(cell_c_chunks{i});
        epsilon_spmd{i} = eps_b(cell_c_chunks{i});
        aW_spmd{i} = aW(cell_c_chunks{i});
        W_spmd{i} = W(cell_c_chunks{i});
        G_spmd{i} = G(cell_c_chunks{i});
    end
    param_HSI.ind = ind;
    clear yb eps_b aW W G

    [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter] = ...
        facetHyperSARA_cw_rd2(y_spmd, epsilon_spmd, ...
        A, At, aW_spmd, G_spmd, W_spmd, param_HSI, Qx, Qy, ncores_data, ...
        wlt_basis, filter_length, nlevel, cell_c_chunks, ch(end), overlap_size, window_type, ...
    fullfile(results_path,warm_start(nChannels)), fullfile(results_path,temp_results_name(nChannels)), flag_nonZeroPrimal, flag_homotopy);

    time_iter_average = mean(end_iter);
    disp(['Average time per iteration: ', num2str(time_iter_average)]);

    %! TO BE CHECKED
    mkdir('results/')
    save(fullfile(results_path, results_name),'-v7.3','xsol', ...
        'param', 'epsilon', 'rel_fval', 'nuclear', 'l21', 'norm_res_out', ...
        'end_iter', 'time_iter_average');
    fitswrite(xsol,fullfile(results_path, strcat('x_fhs_', algo_version, ...
        '_', window_type, ...
        '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc), ...
        '_ind=', num2str(ind), ...
        '_overlap=', num2str(overlap_size), ...
        '_primal=', num2str(flag_nonZeroPrimal), ...
        '_homotopy=', num2str(flag_homotopy), ...
        '.fits')))
    %!----
end
