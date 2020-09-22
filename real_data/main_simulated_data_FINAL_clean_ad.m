function main_simulated_data_FINAL_clean_ad(subcubeIdx,param_global,jobID)
% l1_req: soft thresholding param
% rw : reweights
%% flags
%% read essential params
%:AD:Default params
%%? to be modified ?
% data specific
if ~isfield(param_global, 'visibilityFileName'),  param_global.visibilityFileName = {'5','7'};end
%if ~isfield(param_global, 'telescope'),  param_global.telescope = {'vla'}; end % configurations are specific to VLA
if ~isfield(param_global, 'configurationFileName'),  param_global.configurationFileName = {'C','A'}; end

% image resolution & dimensions
if ~isfield(param_global, 'pixelSize'),   param_global.pixelSize = []; end %in arcsec
if isempty (param_global.pixelSize),  param_global.imageResolution = 'nominal';
else,  param_global.imageResolution = 'user_defined';
end
if ~isfield(param_global, 'Nx'),  param_global.Nx = 2048; end
if ~isfield(param_global, 'Ny'),  param_global.Ny = 2048; end
% facet related
if ~isfield(param_global, 'Qx'),  param_global.Qx =  floor(param_global.Nx/512); end
if ~isfield(param_global, 'Qy'),  param_global.Qy =  floor(param_global.Ny/512); end
%
if ~isfield(param_global, 'Qc2'), param_global.Qc2 = 20; end
if ~isfield(param_global, 'window_type'), param_global.window_type = 'triangular'; end
if ~isfield(param_global, 'parallel_version'), param_global.parallel_version ='spmd4_cst_weighted'; end
% data blocks
if ~isfield(param_global, 'nDataBlk'),  param_global.nDataBlk = []; end
if ~isfield(param_global, 'sizeDataBlk'),  param_global.sizeDataBlk = []; end
if isempty(param_global.sizeDataBlk) && isempty (param_global.nDataBlk)
    param_global.sizeDataBlk = 2e5;
end

% reweighting
if ~isfield(param_global,'currentReweightStep' ), param_global.currentReweightStep = [];  end
% sparsity dict.
if ~isfield(param_global,'wavelet_level' ), param_global.wavelet_level = 4;  end
if ~isfield(param_global,'wavelet_basis' ), param_global.wavelet_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};  end
% w projection
if ~isfield(param_global, 'CEnergyL2') , param_global.CEnergyL2 =1-5e-5; end
if ~isfield(param_global, 'GEnergyL2') , param_global.GEnergyL2 =1-1e-4; end
%  algo related
if ~isfield(param_global,'select_algo' ), param_global.select_algo = 2 ; end
if ~isfield(param_global, 'l1_reg') , param_global.l1_reg =1e-6; end
if ~isfield(param_global, 'rel_obj') , param_global.rel_obj =1e-8; end
if ~isfield(param_global, 'max_iter') , param_global.max_iter = 1e6; end
% flags
if ~isfield(param_global,'flag_wProjection' ), param_global.flag_wProjection = 1 ; end
if ~isfield(param_global,'flag_extract_raw_data' ), param_global.flag_extract_raw_data =0 ; end
if ~isfield(param_global,'flag_compute_Anorm' ),    param_global.flag_compute_Anorm =0 ; end
if ~isfield(param_global,'flag_generate_real_data' ), param_global.flag_generate_real_data = 1 ; end
if ~isfield(param_global,'flag_generate_eps_nnls' ),  param_global.flag_generate_eps_nnls = 1 ; end
if ~isfield(param_global,'flag_save_data' ), param_global.flag_save_data =0  ; end
if ~isfield(param_global,'flag_save_full_operator' ), param_global.flag_save_full_operator =0  ; end
if ~isfield(param_global,'flag_free_memory'), param_global.flag_free_memory =1 ; end
if ~isfield(param_global,'flag_solve_minimization'), param_global.flag_solve_minimization =1 ; end
if ~isfield(param_global,'flag_load_temp_result' ),   param_global.flag_load_temp_result =0 ; end
if ~isfield(param_global,'flag_plotFigures' ), param_global.flag_plotFigures = 0;  end
if ~isfield(param_global,'flag_verbose' ), param_global.flag_verbose = 0;  end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% get flags
doExtractRawData    = param_global.flag_extract_raw_data;
doGenerateRealData  = param_global.flag_generate_real_data;
doWProjection = param_global.flag_wProjection;

%:AD:  we need to add the theortical option for setting the l2 bounds (chi_squared distribution)
doGenerateEpsNNLS   = param_global.flag_generate_eps_nnls;
%
doComputeAnorm       = param_global.flag_compute_Anorm;
doSolveMinimization  = param_global.flag_solve_minimization;
%
doLoadTempResult    = param_global.flag_load_temp_result;
doSaveFullOperator  = param_global.flag_save_full_operator;
doSaveData   = param_global.flag_save_data;
doFreeMemory = param_global.flag_free_memory;
doPlotFigures = param_global.flag_plotFigures;
verbose =  param_global.flag_verbose;
% data
vlaConfigFileName = param_global.configurationFileName;
visibilityFileName = param_global.visibilityFileName;

% image dimensions & resolution
Nx = param_global.Nx;
Ny = param_global.Ny;
pixelSize  = param_global.pixelSize;
imageResolutionStr = param_global.imageResolution;
switch imageResolutionStr
    case 'nominal'
        fprintf('\nWARNING: No pixelsize provided by user --> adopting 1.5x instrumental resolution.\n')
    otherwise
        fprintf('\nINFO: Pixelsize provided by user (in arcsec by default): %f asec.\n',pixelSize);
end

% meas. op.: data blocks
% data blocks
nDataBlk= param_global.nDataBlk ;
szDataBlk= param_global.sizeDataBlk ;

% meas. op.: w projection
wproj.CEnergyL2 = param_global.CEnergyL2 ;
wproj.GEnergyL2 = param_global.GEnergyL2 ;

% wavelet related
wvlt.nlevel = param_global.wavelet_level ; % wavelet level
wvlt.basis  = param_global.wavelet_basis; % wavelet basis to be used, always put self in last position if used
wvlt.FilterLen = [2*(1:8)'; 0]; % length of the filters (0 corresponding to the 'self' basis)

% facet related
Qx = param_global.Qx;
Qy = param_global.Qy;
Qc2 = param_global.Qc2;
window_type = param_global.window_type;
parallel_version = param_global.parallel_version;
if strcmp(parallel_version, 'spmd4_cst_weighted') || strcmp(parallel_version, 'spmd4_cst')
    %:AD: this var should be renamed
    d = 512; %(power(2, nlevel)-1)*(2*8 - 2); % assuming db8 largest wavelet filter
end

% reweighting
currentReweightStep = param_global.currentReweightStep;
% algo related
param_algo.select_algo  = param_global.select_algo;
param_algo.l1_reg = param_global.l1_reg;
param_algo.rel_obj  = param_global.rel_obj; % stopping criterion
param_algo.max_iter = param_global.max_iter ; % max number of iterations
%% setting paths
%:AD: to be modified
pathProject = pwd;
irtLibrary = [pathProject,'/nufft/'];
addpath(irtLibrary)
addpath(genpath([pathProject,'/real_data/']));
addpath(genpath([pathProject,'/hypersara-clean/lib/']));
addpath(genpath([pathProject,'/lib/']));
addpath(genpath([pathProject,'/sdwt2/']));
addpath(genpath([pathProject,'/src/']))
addpath(genpath([pathProject,'/data/']))
pathData =  [pathProject,'/data/'];
pathResults =  [pathProject,'/results/'];
mkdir(pathResults)

%% Default NUFFT params (invoked in Fessler's code)
nufft.ox = 2; % oversampling factors for nufft
nufft.oy = 2; % oversampling factors for nufft
%:AD: can we limit the kernel to 7x7 using Kaiser Bessel only,for a smaller G matrix
%%? no minmax approachin c++ anyway ?
nufft.Kx = 8; % number of neighbours for nufft
nufft.Ky = 8; % number of neighbours for nufft

nFourier  = prod([nufft.oy*Ny nufft.ox*Nx]) ;

%% Generate or Load real data
if doGenerateRealData
    
    %% preconditioning param
    param_precond.N = Nx*Ny; % number of pixels in the image
    param_precond.Nox = nufft.ox*Nx; % number of pixels in the image
    param_precond.Noy = nufft.oy*Ny; % number of pixels in the image
    param_precond.gen_uniform_weight_matrix = 1; %set weighting type
    param_precond.uniform_weight_sub_pixels = 1;
    
    %% block structure param
    %:AD: extra flags to be added for partionning
    
    %     regenerate_block_structure = 1;
    param_block_structure.use_density_partitioning = 0;
    param_block_structure.density_partitioning_no = 1;
    %
    param_block_structure.use_uniform_partitioning = 0;
    param_block_structure.uniform_partitioning_no = 4;
    %
    param_block_structure.use_equal_partitioning = 0;
    param_block_structure.equal_partitioning_no = 1;
    %
    param_block_structure.use_manual_frequency_partitioning = 0;
    % sparam.fpartition = [pi]; % partition (symetrically) of the data to nodes (frequency ranges)
    % sparam.fpartition = [0, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
    % sparam.fpartition = [-0.25*pi, 0, 0.25*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
    % sparam.fpartition = [-64/256*pi, 0, 64/256*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
    param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi/4), 0, icdf('norm', 0.75, 0, pi/4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)
    % sparam.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
    % sparam.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
    
    param_block_structure.use_manual_partitioning = 1;
    
    %% load real_data;
    if doExtractRawData % MODIFY READING HERE (MATFILE)
        [y, uw1, vw1,ww, nWw, f, time_, pos] = util_load_real_vis_data_ad(visibilityFileName,vlaConfigFileName, subcubeIdx,pathData);
        [uw, vw,pixelSize] = scale_uv_ad(visibilityFileName,vlaConfigFileName, uw1, vw1, pixelSize,imageResolutionStr,pathData);
        %
        save([pathData '/CYG_data_raw_ind=' num2str(subcubeIdx) '.mat'],'-v7.3', 'y', 'uw', 'vw','ww', 'nWw', 'f', 'time_', 'pos', 'pixelSize');
    else
        %:AD: modified path
        load([pathData '/CYG_data_raw_ind=' num2str(subcubeIdx) '.mat'],'y', 'uw', 'vw', 'ww', 'nWw', 'time_', 'pos','pixelSize');
        fprintf('\nINFO: reading data from a saved mat file..\nINFO: uv coor. are scaled in [-pi pi].')
        fprintf('\nINFO: pixelsize: %f asec\n',pixelSize)
    end
    
    %% FoV info 
    FoVx = pixelSize * Nx *pi/180/3600;
    FoVy = pixelSize * Ny *pi/180/3600;
    uvGridSizex   = 1/(nufft.ox*FoVx);
    uvGridSizey   = 1/(nufft.oy*FoVy);
    minGridSize   = min(uvGridSizey,uvGridSizex); %smallest size of the gridcell in the spatial Fourier domain
    
    %% get operators
    nChannels  = size(y,2);  
    channels = [1 : nChannels]; % 

    if doPlotFigures
        figure(1), %AD: plot uv coverage for all channels
        for iCh = 1:nChannels
            scatter(uw{iCh},vw{iCh},'.');hold on;
        end
    end
  
    G  = cell(nChannels,1);
    aW = cell(nChannels,1);
    W  = cell(nChannels,1);
    yb = cell(nChannels,1);
    
    for iCh = channels
       
        nMeasPerCh = length(uw{iCh});

        %% set the blocks structure
        if param_block_structure.use_manual_partitioning == 1
            if ~isempty(nDataBlk),   szDataBlk = floor(nMeasPerCh/nDataBlk);
            else, szDataBlk = szDataBlk*(nMeasPerCh>=2*szDataBlk) + nMeasPerCh * (nMeasPerCh< 2*szDataBlk);
            end
            param_block.size = szDataBlk;
            param_block.snapshot = 0;
            param_block.pos = pos{iCh};
            out_block = util_time_based_block_sp_ar(uw{iCh},time_{iCh},param_block);
            partition = out_block.partition;
            param_block_structure.partition = partition;
        end
        %% compute weights
        aWw = util_gen_preconditioning_matrix(uw{iCh}, vw{iCh}, param_precond);
        % set the blocks structure
        [u, v, ~, uvidx, aW{iCh}, nW] = util_gen_block_structure(uw{iCh}, vw{iCh}, aWw, nWw{iCh}, param_block_structure); % CHANGE DEF. OF OPERATOR
        uw{iCh} = [];
        vw{iCh} = [];
        nWw{iCh} = [];
        %% measurement operator initialization
        fprintf('Initializing the NUFFT operator\n\n');
        tstart = tic;
        %
        [A, At, G{iCh}, W{iCh}] = op_p_nufft([v u], [Ny Nx], [nufft.Ky nufft.Kx], [nufft.oy*Ny nufft.ox*Nx], [Ny/2 Nx/2], nW);
        %% restructure data
        nBlk    = length(u);
        yb{iCh} = cell(nBlk,1);
        for jBlk = 1 : nBlk,   yb{iCh}{jBlk} = y{iCh}(uvidx{jBlk});
        end
        
        %% w-correction
        %:AD:  check if w correction is necessary
        if ~doWProjection
            effBandwidthWterm = max(abs(max(FoVy, FoVx).*ww{iCh}(:)));
            if effBandwidthWterm > 3*minGridSize
                fprintf('\nINFO: W correction will be performed\n')
                doWProjection = 1;
            else
                fprintf('\nINFO: NO W correction\n')
                doWProjection = 0;
            end
        end
        
        if  doWProjection
            w   = cell(nBlk,1);
            for jBlk = 1 : nBlk,   w{jBlk} =  ww{iCh}(uvidx{jBlk});
            end
            ww{iCh} =[];
            param_wterm.FoV =[FoVy; FoVx];
            param_wterm.ox = [nufft.oy ;nufft.ox];
            param_wterm.gImDims = [Ny; Nx];
            
            for jBlk = 1:nBlk
                %:AD: G back to sparse matrix
                GBis = sparse(size(G{iCh}{jBlk},1), nFourier);
                GBis(:,W{iCh}{jBlk}) =  G{iCh}{jBlk};
                %:AD: build the w proj. op.
                GBis =  getWprojGmatrix(GBis, w{jBlk},param_wterm,wproj.CEnergyL2,wproj.GEnergyL2);
                %:AD: restructure G matrix
                W{iCh}{jBlk} = false(nFourier, 1);
                W{iCh}{jBlk} = any(GBis, 1).';
                G{iCh}{jBlk}= GBis(:, W{iCh}{jBlk});
                clear G_dummy;
            end
        end
        
    end
    
    %% Save data
    if doSaveData
        save([pathData 'CYG_data.mat'],'-v7.3', 'G', 'W', 'aW', 'yb');
        save([pathData 'CYG_y.mat'],'-v7.3', 'y');
    end
    %
    if doSaveFullOperator && exist('Gw','var')
        save([pathData 'CYG_Gw.mat'],'-v7.3', 'Gw');
    end
    
    
    %% Free memory
    if doFreeMemory
        clear y u v uv_mat uv_mat1 uw vw ww ant1 ant2 aWw nW nWw out_block;
    end
    % % %     extract_real_data;
else
    load([pathData 'CYG_data.mat'],'G', 'W', 'aW', 'yb');
    nChannels =  size(yb,2);
    channels = [1 : nChannels];

    [A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [nufft.Ky nufft.Kx], [nufft.oy*Ny nufft.ox*Nx], [Ny/2 Nx/2]);
end

%% Generate facets
% spectral faceting (interlaved sampling)
%id = interleaved_facets(Qc, c);

% spatial faceting
Q = Qx*Qy;
rg_y = domain_decomposition(Qy, Ny);
rg_x = domain_decomposition(Qx, Nx);

% create starting index of the spatial facets (I) and associated dimensions
% (dims). Beware: I starts at (0, 0)
I = zeros(Q, 2);
dims = zeros(Q, 2);
for qx = 1:Qx
    for qy = 1:Qy
        q = (qx-1)*Qy+qy;
        I(q, :) = [rg_y(qy, 1)-1, rg_x(qx, 1)-1];
        dims(q, :) = [rg_y(qy,2)-rg_y(qy,1)+1, rg_x(qx,2)-rg_x(qx,1)+1];
    end
end

%% Compute map estimator
if doSolveMinimization
    % % % % %     Solver_simulated_data_FINAL_clean_ad
    
    %% Generate initial epsilons by performing imaging with NNLS on each data block separately
    if doGenerateEpsNNLS
        % param_nnls.im = im; % original image, used to compute the SNR
        param_nnls.verbose = verbose; % print log or not
        param_nnls.rel_obj = 1e-3; % stopping criterion
        param_nnls.max_iter = 1000; % max number of iterations
        param_nnls.sol_steps = [inf]; % saves images at the given iterations
        param_nnls.beta = 1;
        param_nnls.nFourier = nFourier;
        
        % solve nnls per block
        % util_create_pool(min(channels,24));
        eps_b = cell(nChannels,1);
        for iCh = channels
            nBlk = length(G{iCh});
            eps_b{iCh} = cell(nBlk,1);
            %:AD
            fprintf('\nSolving NNLS for channel  %i\n',iCh)
            G_band = [];
            partitionBis = 1;
            for jBlk = 1 : nBlk
                GBis = sparse(size(G{iCh}{jBlk},1), nFourier);
                GBis(:,W{iCh}{jBlk}) =  G{iCh}{jBlk};
                G_band = sparse([G_band ;GBis]);
                partitionBis = [partitionBis size(G{iCh}{jBlk},1)];
            end
            [~,resNNLS,~] = fb_nnls_blocks_ad(vertcat(yb{iCh}{:}), A, At, G_band, ':', param_nnls);
            for jBlk = 1 : nBlk
                eps_b{iCh}{jBlk} =  norm(resNNLS(sum(partitionBis(1:jBlk)):sum(partitionBis(1:jBlk+1))-1));
            end
            clear  resNNLS  GBis  G_band partitionBis;
            
            %%:AD: previous code to compute l2 per blk
            %   for j = 1 : length(G{i})
            %        % printf('solving for band %i\n\n',i)
            %        [~,eps_b{i}{j}] = fb_nnls_blocks(vertcat(yb{i}{j}, A, At, G{i}{j}, W{i}{j}, param_nnls);
            %   end
        end
        
        save([pathResults 'eps_ind=' num2str(subcubeIdx) '.mat'],'-v7.3', 'eps_b');
        
    else  
        load([pathResults 'eps_ind=' num2str(subcubeIdx) '.mat'],'eps_b');
    end
    
    %% compute op. norm
    if doComputeAnorm
        %Compute full measurement operator spectral norm
        F = afclean( @(x) HS_forward_operator_precond_G(x, G, W, A, aW));
        Ft = afclean( @(y) HS_adjoint_operator_precond_G(y, G, W, At, aW, Ny, Nx));
        Anorm = getSpectralNormPowerMethod_ad(F, Ft, [Ny Nx length(channels)],1e-8,1000,verbose);
        save([pathResults 'Anorm_ind=' num2str(subcubeIdx) '.mat'],'-v7.3', 'Anorm');
    else,   load([pathResults 'Anorm_ind=' num2str(subcubeIdx) '.mat'], 'Anorm');
    end
    
    
    %% Splitting operator FULL
    % Sp = @(x) Split_forward_operator(x,I,dims,Q);
    % Spt = @(x) Split_adjoint_operator(x,I,dims,Q,Ny,Nx,length(ch));
    % Sp_norm = pow_method_op(Sp,Spt,[Ny Nx length(ch)]); % [P.-A.] in theory, Sp_pnorm = Pnorm (no need to compute both)
    
    %% HSI parameter structure sent to the  HSI algorithm
    % user defined
    param_HSI.gamma = param_algo.l1_reg;  %convergence parameter L1 (soft th parameter)
    param_HSI.rel_obj = param_algo.rel_obj; % stopping criterion
    param_HSI.max_iter = param_algo.max_iter; % max number of iterations
    %
    param_HSI.verbose = 2; % print log or not
    param_HSI.nu0 = 1; % bound on the norm of the Identity operator
    param_HSI.nu1 = 1; % bound on the norm of the operator Psi
    param_HSI.gamma0 = 1;
    %
    param_HSI.use_adapt_eps = 1; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
    param_HSI.adapt_eps_start = 300; % minimum num of iter before stating adjustment
    param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
    param_HSI.adapt_eps_tol_out = 1.001; % tolerance outside the l2 ball
    param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
    param_HSI.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
    param_HSI.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data
    
    param_HSI.reweight_alpha = 1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
    param_HSI.reweight_alpha_ff = 0.8; 0.9;
    param_HSI.total_reweights = 70; -1; % -1 if you don't want reweighting
    param_HSI.reweight_abs_of_max = 1; 0.005; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised
    
    param_HSI.use_reweight_steps = 1; % reweighting by fixed steps
    param_HSI.reweight_step_size = 300; % reweighting step size
    param_HSI.reweight_steps = [5000: param_HSI.reweight_step_size :10000];
    param_HSI.step_flag = 1;
    
    param_HSI.use_reweight_eps = 0; % reweighting w.r.t the relative change of the solution
    param_HSI.reweight_max_reweight_itr = param_HSI.max_iter - param_HSI.reweight_step_size;
    param_HSI.reweight_rel_obj = 1e-4; % criterion for performing reweighting
    param_HSI.reweight_min_steps_rel_obj = 300; % min num of iter between reweights
    
    param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid
    
    
    %% HyperSARA-sdwt2 (split L21 + nuclear norms + wavelets)
    
    if param_algo.select_algo ==2
        param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
        disp('Split L21 + Nuclear + wavelets')
        % spectral tesselation (non-overlapping)
        rg_c = domain_decomposition(Qc2, channels(end));
        cell_c_chunks = cell(Qc2, 1);
        y_spmd = cell(Qc2, 1);
        epsilon_spmd = cell(Qc2, 1);
        aW_spmd = cell(Qc2, 1);
        W_spmd = cell(Qc2, 1);
        G_spmd = cell(Qc2, 1);
        
        for iCh = 1:Qc2
            cell_c_chunks{iCh} = rg_c(iCh, 1):rg_c(iCh, 2);
            y_spmd{iCh} = yb(cell_c_chunks{iCh});
            epsilon_spmd{iCh} = eps_b(cell_c_chunks{iCh});
            aW_spmd{iCh} = aW(cell_c_chunks{iCh});
            W_spmd{iCh} = W(cell_c_chunks{iCh});
            G_spmd{iCh} = G(cell_c_chunks{iCh});
        end
        
        clear yb eps_b aW W G;
        
        if  currentReweightStep >= 0
            param_HSI.ind = subcubeIdx;
            load([pathResults '/result_HyperSARA_spmd4_cst_weighted_rd_' ...
                num2str(param_HSI.ind) '_' num2str(param_HSI.gamma) '_' num2str(currentReweightStep) '.mat'],'epsilon','param');
            epsilon_spmd = epsilon;
            param_HSI.init_xsol = param.init_xsol;
            param_HSI.init_g  = param.init_g;
            param_HSI.init_v0 = param.init_v0;
            param_HSI.init_v1 = param.init_v1;
            param_HSI.init_weights0 = param.init_weights0;
            param_HSI.init_weights1 = param.init_weights1;
            param_HSI.init_v2 = param.init_v2;
            param_HSI.init_proj = param.init_proj;
            param_HSI.init_t_block = param.init_t_block;
            param_HSI.init_t_start = param.init_t_start+1;
            param_HSI.reweight_alpha = param.reweight_alpha;
            param_HSI.init_reweight_step_count = param.init_reweight_step_count;
            param_HSI.init_reweight_last_iter_step = param.init_reweight_last_iter_step;
            param_HSI.reweight_steps = (param_HSI.init_t_start+param_HSI.reweight_step_size: param_HSI.reweight_step_size :param_HSI.max_iter+(2*param_HSI.reweight_step_size));
            param_HSI.step_flag = 0;
        end
        
        switch parallel_version
            case 'spmd4_cst_weighted'% same as spmd_sct, weight correction (apodization window in this case)
                [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
                    pdfb_LRJS_precond_NL21_sdwt2_spmd4_cst_overlap_weighted_rd(y_spmd, epsilon_spmd, ...
                    A, At, aW_spmd, G_spmd, W_spmd, param_HSI, Qx, Qy, Qc2, ...
                    wvlt.basis, wvlt.FilterLen, wvlt.nlevel, cell_c_chunks, channels(end), d, window_type);
                
            otherwise
                error('Unknown parallelisation option.')
        end
        
        averageTimePerItr = mean(end_iter)
        
        save([pathResults,'/results_hyperSARA_', parallel_version, '_ind=',...
            num2str(subcubeIdx), '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', ...
            num2str(Qc2), '_gamma=', num2str(gamma),'.mat'],'-v7.3',...
            'xsol', 'param', 'epsilon', 't', 'rel_fval', 'nuclear', 'l21', 'norm_res', 'res','end_iter');
        
    end
    
    
    
    
    
end

end
