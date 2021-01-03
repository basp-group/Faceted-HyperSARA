clear all;
close all;
clc;
format compact;

addpath lib/generate_data/
addpath lib/operators/
addpath lib/measurement-operator/nufft/
addpath lib/utils/
addpath lib/sdwt2/
addpath data/
addpath src/
addpath src/spmd
addpath src/spmd/weighted % standard no_overlap weighted: to be added depending on the configuration
addpath src/spmd/standard
addpath src/spmd/no_overlap
% addpath lib/irt
% 
% try
%     run('lib/irt/setup.m');
% end

%% /!\ REVISE DEFAULT PARAMETERS /!\

ind = 1; % index from "spectral" facet
img_size = 512;
Qc2 = 1;
% T = 1500; % to be set depending on the value of p
% hrs = 6;

flag_algo = 2;
img_name = 'W28_512';
Nx = 512;
Ny = Nx;
p = 0.3; % percentage
kernel = 'minmax:tuned'; % 'kaiser', 'minmax:tuned'
nchannels = 60; % nchannels
generate_cube = true;
generate_coverage = true;
generate_visibilities = true;
compute_operator_norm = true;

generate_undersampled_cube = true; % Default 15 channels cube with line emissions
solve_minimization = true;
generate_eps_nnls = 0;

% parallel_version = % options: 
Qx = 2;
Qy = 1;
Qc = 1;
input_snr = 40; % input SNR (in dB)
% num_chunk
cube_path = strcat('data/', img_name, '_L', num2str(nchannels), '.fits');
coverage_path = 'data/uv_coverage_p=1.mat';
window_type = 'triangular';
parallel_version = 'cst_weighted';
save_data = 1; % keep to 1 by default

data_path = 'data/'; % path to data, image cube, coverage 
results_path = 'results/';
path_reference_cube = fullfile(data_path, strcat(img_name, '.fits'));
data_name_function = @(nchannels) strcat('y_N=',num2str(Nx),'_L=',num2str(nchannels),'_p=',num2str(p),'_snr=', num2str(input_snr),'.mat');
results_name_function = @(nchannels) strcat('results_N=',num2str(Nx),'_L=',num2str(nchannels),'_p=',num2str(p),'_snr=', num2str(input_snr),'.mat');

% max_iter
% max_reweight

% select algorithm parameters
bool_weights = true; % for the spmd4_new version (50% overlap version)

if strcmp(parallel_version, 'cst_weighted') || strcmp(parallel_version, 'cst')
    nlevel = 4;
    d = (power(2, nlevel)-1)*(2*8 - 2); % assuming db8 largest wavelet filter
end

%% Default config parameters (can be kept as is)
% gridding parameters
N = Nx * Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

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
param_block_structure.use_equal_partitioning = 1;
param_block_structure.equal_partitioning_no = 1;
param_block_structure.use_manual_frequency_partitioning = 0;
% sparam.fpartition = [pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [0, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.25*pi, 0, 0.25*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-64/256*pi, 0, 64/256*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_block_structure.fpartition = [icdf('norm', 0.25, 0, pi/4), 0, icdf('norm', 0.75, 0, pi/4), pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
param_block_structure.use_manual_partitioning = 0;

% setting random number generator (if data are generated)
rng(1234)

mkdir(data_path)
mkdir(results_path)

%% Generate/load ground-truth image cube
f = linspace(1,2,nchannels);
if generate_cube
    % frequency bandwidth from 1 to 2 GHz
    emission_lines = 0; % insert emission lines on top of the continuous spectra
    [x0,X0] = Generate_cube(path_reference_cube,f,emission_lines);
    [Ny, Nx, nchannels] = size(x0);
    if generate_undersampled_cube
        % undersampling factor for the channels
        unds = 4; % take 1/unds images
        [x0,X0,f,nchannels] = Generate_undersampled_cube(x0,f,Ny,Nx,nchannels,unds);
    end
    fitswrite(X0.', cube_path);
    fitsdisp(cube_path)
else
    X0 = fitsread(cube_path).';
    x0 = reshape(X0, [Ny, Nx, nchannels]);
end

data_name = data_name_function(nchannels);
results_name = results_name_function(nchannels);

%% Generate facets
% spectral faceting (interleaved sampling)
id = interleaved_facets(Qc, nchannels);

% create the complete tessellation (spectral)
if ind > 0
    x0 = x0(:,:,id{ind});
    c = size(x0,3);
    channels = 1:c;
    f = f(id{ind});
    X0 = reshape(x0,Nx*Ny,c);
end

%% /!\ SECTION TO BE REVISED /!\
if flag_algo == 0 % L11
    param_HSI.num_workers = 34;
elseif flag_algo == 1 % HyperSARA
    param_HSI.num_workers = 1;
elseif flag_algo == 2 % Faceted HyperSARA
    param_HSI.num_workers = Qx*Qy*2+1;  %%%%%%%%%%% TO BE SET BY P.-A.
end
param_HSI.num_workers

%% Generate/load uv-coverage, setup measurement operator
% generating u-v coverage
if generate_coverage
    cov_type = 'vlaa';
    dl = 1.1;
    hrs = 5;
    na = 27; % for vlaa
    M = na*(na-1)/2;
    % Fixing Mt = 0.5 N, take T = 0.5 N / M : na = 27 for vla
    T = floor(p*(Nx*Ny)/M); % should be > 1
    [u, v, ~] = generate_uv_coverage(T, hrs, dl, cov_type);
    % -> in the c++ code, do u.col(f) = u*freq[f]/freq[0]
    % with freq = linspace(1, 2, L)
    % om = [v(:), u(:)];
    u = u(:)*f(1)/f(end);
    v = v(:)*f(1)/f(end);
    fitswrite([u, v, ones(numel(u), 1)], coverage_path)
    fitsdisp(coverage_path);
else
    uvw = fitsread(coverage_path);
    u = uvw(:, 1);
    v = uvw(:, 2);
    % u = u*f(1)/f(end);
    % v = v*f(1)/f(end);
    disp('Coverage loaded successfully')
    clear uvw
end

% setup measurement operator
for i = 1:nchannels
    uw = (f(i)/f(1)) * u;
    vw = (f(i)/f(1)) * v;
    
    % compute uniform weights (sampling density) for the preconditioning
    aWw = util_gen_preconditioning_matrix(uw, vw, param_precond);
    
    % set the weighting matrix, these are the natural weights for real data
    nWw = ones(length(uw), 1);
    
    % set the blocks structure
    [u1, v1, ~, ~, aW{i}, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    [A, At, G{i}, W{i}] = op_p_nufft([v1 u1], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
end

%% Free memory
clear u v u1 v1 uw vw aWw nW nWw param_block_structure param_precond;

%% Generate/load visibilities
if generate_visibilities
    param_l2_ball.type = 'sigma';
    param_l2_ball.sigma_ball = 2;
    [y0, y, Nm, sigma_noise] = util_gen_measurements(x0, G, W, A, input_snr);
    [epsilon,epsilons] = util_gen_data_fidelity_bounds(y, Nm, param_l2_ball, sigma_noise); 

    if save_data
        save(fullfile(data_path,data_name), '-v7.3', 'y0', 'y', 'epsilon', 'epsilons', 'sigma_noise');
    end
else
    load(fullfile(data_path,data_name), 'y', 'epsilon', 'epsilons', 'sigma_noise');
    disp('Data loaded successfully')
end

clear y0 Nm;

%% Compute MAP estimator
if solve_minimization
    if compute_operator_norm
        if flag_algo == 0 % SARA solver
            % Compute measurement operator spectral norm for each channel individually
            Anorm = zeros(nchannels, 1);
            for i = channels
                %Anorm_ch(i) = pow_method_op(@(x) sqrt(cell2mat(aW{i})) .* (Gw{i}*A(x)), @(x) real(At(Gw{i}' * (sqrt(cell2mat(aW{i})) .* x))), [Ny Nx 1]);
                F = afclean( @(x) HS_forward_operator_precond_G(x, G(i), W(i), A, aW(i)));
                Ft = afclean( @(y) HS_adjoint_operator_precond_G(y, G(i), W(i), At, aW(i), Ny, Nx));
                Anorm(i) = pow_method_op(F, Ft, [Ny Nx 1]);
                save(fullfile(results_path,strcat('Anorm_SARA_N=',num2str(Nx), ...
            '_L=',num2str(nchannels),'_p=',num2str(p),'_snr=', num2str(input_snr), '.mat')),'-v7.3', 'Anorm');
            end
        else
            % Compute full measurement operator spectral norm
            F = afclean( @(x) HS_forward_operator_precond_G(x, G, W, A, aW));
            Ft = afclean( @(y) HS_adjoint_operator_precond_G(y, G, W, At, aW, Ny, Nx));
            Anorm = pow_method_op(F, Ft, [Ny Nx length(channels)]); 
            save(fullfile(results_path,strcat('Anorm_N=',num2str(Nx), ...
            '_L=',num2str(nchannels),'_p=',num2str(p),'_snr=', num2str(input_snr), '.mat')),'-v7.3', 'Anorm');
        end
    else
        if flag_algo == 0
            load(fullfile(results_path,strcat('Anorm_SARA_N=',num2str(Nx),'_L=', ...
                num2str(nchannels),'_p=',num2str(p),'_snr=', num2str(input_snr), '.mat')));
        else
            load(fullfile(results_path,strcat('Anorm_N=',num2str(Nx),'_L=', ...
                num2str(nchannels),'_p=',num2str(p),'_snr=', num2str(input_snr), '.mat')));
        end
    end
        
    
    %% Generate initial epsilons by performing imaging with NNLS on each data block separately
    if generate_eps_nnls
        % param_nnls.im = im; % original image, used to compute the SNR
        param_nnls.verbose = 2; % print log or not
        param_nnls.rel_obj = 1e-5; % stopping criterion
        param_nnls.max_iter = 1000; % max number of iterations
        param_nnls.sol_steps = [inf]; % saves images at the given iterations
        param_nnls.beta = 1;
        % solve nnls per block
        for i = channels
            eps_b{i} = cell(length(G{i}),1);
            for j = 1 : length(G{i})
                % printf('solving for band %i\n\n',i)
                [~,eps_b{i}{j}] = fb_nnls_blocks(yb{i}{j}, A, At, G{i}{j}, W{i}{j}, param_nnls);
            end
        end
        if save_data
            mkdir('data/')
            save('data/eps.mat','-v7.3', 'eps_b');
        end
    elseif load_data
        load('data/eps.mat');
    end
    
    %% sparsity operator definition
    nlevel = 4; % wavelet level
    wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used, always put self in last position if used
    L = [2*(1:8)'; 0]; % length of the filters (0 corresponding to the 'self' basis)
    
    if flag_algo < 2
        
        [Psi1, Psit1] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
        P = length(Psi1);
        
        for k = 1 : P
            f = '@(x_wave) HS_forward_sparsity(x_wave,Psi1{';
            f = sprintf('%s%i},Ny,Nx);', f,k);
            Psi{k} = eval(f);
            
            b(k) = size(Psit1{k}(zeros(Ny,Nx,1)),1);
            
            ft = strcat('@(x) HS_adjoint_sparsity(x,Psit1{', num2str(k), '},b(', num2str(k), '));');
            Psit{k} = eval(ft);
        end
        
        %% Full sparsity operator
        [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
        bb = size(Psitw(zeros(Ny,Nx,1)),1);
        
        Psi_full = @(x_wave) HS_forward_sparsity(x_wave,Psiw,Ny,Nx);
        Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,bb);
    end
    
    %% HSI parameter structure sent to the  HSI algorithm
    param_HSI.verbose = 2; % print log or not
    param_HSI.nu0 = 1; % bound on the norm of the Identity operator
    param_HSI.nu1 = 1; % bound on the norm of the operator Psi
    param_HSI.gamma0 = 1;
    param_HSI.gamma = 1e-2;  %convergence parameter L1 (soft th parameter)
    param_HSI.rel_obj = 1e-5; % stopping criterion
    param_HSI.max_iter = 10000; % max number of iterations
    
    param_HSI.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
    param_HSI.adapt_eps_start = 200; % minimum num of iter before stating adjustment
    param_HSI.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
    param_HSI.adapt_eps_tol_out = 1.001; % tolerance outside the l2 ball
    param_HSI.adapt_eps_steps = 100; % min num of iter between consecutive updates
    param_HSI.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
    param_HSI.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data
    
    param_HSI.reweight_alpha = 1; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
    param_HSI.reweight_alpha_ff = 0.9;
    param_HSI.total_reweights = 30; % -1 if you don't want reweighting
    param_HSI.reweight_abs_of_max = 1; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised
    
    param_HSI.use_reweight_steps = 1; % reweighting by fixed steps
    param_HSI.reweight_step_size = 300; % reweighting step size
    param_HSI.reweight_steps = [5000: param_HSI.reweight_step_size :10000];
    param_HSI.step_flag = 1;
    
    param_HSI.use_reweight_eps = 0; % reweighting w.r.t the relative change of the solution
    param_HSI.reweight_max_reweight_itr = param_HSI.max_iter - param_HSI.reweight_step_size;
    param_HSI.reweight_rel_obj = 1e-4; 5e-4; % criterion for performing reweighting
    param_HSI.reweight_min_steps_rel_obj = 300; % min num of iter between reweights
    
    param_HSI.elipse_proj_max_iter = 20; % max num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_min_iter = 1; % min num of iter for the FB algo that implements the preconditioned projection onto the l2 ball
    param_HSI.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid
    
    %% SARA
    if flag_algo == 0
        param_HSI2 = param_HSI;
        param_HSI2.use_reweight_steps = 0;
        param_HSI2.total_reweights = 30;
        util_create_pool(param_HSI.num_workers);
        maxNumCompThreads(param_HSI.num_workers);
        parfor i = 1 : length(channels)
            [xsol(:,:,i),v1,v2,g,weights1,proj,t_block,reweight_alpha,epsilon,iterh,rel_fval,l11,norm_res,res(:,:,i),end_iter{i}] = ...
                pdfb_L11_Adapt_blocks_rw_precond_new_sim(y(i),epsilons(i), A, At, aW(i), G(i), W(i), Psi_full, Psit_full, param_HSI2, X0(:,i), Anorm(i));
        end
        
        c = size(xsol,3);
        sol = reshape(xsol(:),numel(xsol(:))/c,c);
        snr_x = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)))
        psnrh = zeros(c,1);
        for i = channels
            psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
            Time_iter = mean(end_iter{i});
        end
        snr_x_average = mean(psnrh)
        time_iter_average = mean(Time_iter)
        
        mkdir('results')
        save(['results/result_L11_' num2str(T) '_' num2str(hrs) '.mat'],'-v7.3','xsol', 'sol', 'X0', 'SNR', 'SNR_average', 'res','end_iter', 'Time_iter_average');
        fitswrite(xsol,'results/xsol_L11.fits')
        fitswrite(x0,'results/x0.fits')
    end
    
    %% HyperSARA
    if flag_algo == 1
        param_HSI.nu2 = Anorm; % bound on the norm of the operator A*G
        param_HSI.reweight_alpha_ff = 0.9;
        param_HSI.reweight_abs_of_max = 0.005;
        param_HSI.use_reweight_steps = 0;
        param_HSI.total_reweights = 30;
        param_HSI.use_reweight_eps = 0;
        
        disp('-----------------------------------------')
        disp('HyperSARA')
        disp('-----------------------------------------')
        
        %[xsol,v0,v1,v2,g,weights0,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
        %    pdfb_LRJS_Adapt_blocks_rwNL21_precond_new_sim(y, epsilons, A, At, aW, G, W, Psi_full, Psit_full, param_HSI2, X0);
        
        % run HyperSARA as a special case of Faceted HyperSARA (sub-optimal)
        rg_c = domain_decomposition(Qc2, channels(end));
        cell_c_chunks = cell(Qc2, 1);
        y_spmd = cell(Qc2, 1);
        epsilon_spmd = cell(Qc2, 1);
        aW_spmd = cell(Qc2, 1);
        W_spmd = cell(Qc2, 1);
        G_spmd = cell(Qc2, 1);
        
        for i = 1:Qc2
            cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
            y_spmd{i} = y(cell_c_chunks{i});
            epsilon_spmd{i} = epsilons(cell_c_chunks{i});
            aW_spmd{i} = aW(cell_c_chunks{i});
            W_spmd{i} = W(cell_c_chunks{i});
            G_spmd{i} = G(cell_c_chunks{i});
        end
        
        clear y epsilon aW W G
        
        [xsol,param_HSI,epsilon,t,rel_fval,nuclear,l21,norm_res,res,end_iter] = ...
            pdfb_LRJS_precond_NL21_sdwt2_spmd_serial_SARA(y_spmd, epsilon_spmd, A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qc2, wlt_basis, nlevel, cell_c_chunks, channels(end));
        % [12/06/2019] ok
        
        c = size(xsol,3);
        sol = reshape(xsol(:),numel(xsol(:))/c,c);
        snr_x = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)))
        psnrh = zeros(c,1);
        for i = channels
            psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
        end
        snr_x_average = mean(psnrh)
        time_iter_average = mean(end_iter)
        
        mkdir('results')
        save(['results/result_HyperSARA_' num2str(T) '_' num2str(hrs) '.mat'],'-v7.3','xsol', 'sol', 'X0', 'SNR', 'SNR_average', 'res','end_iter');
        fitswrite(xsol,'results/xsol_HyperSARA.fits')
        fitswrite(x0,'results/x0.fits')
    end    
    
    %% HyperSARA-sdwt2 (split L21 + nuclear norms + wavelets) 
    if flag_algo==2
        
        param_HSI.nu2 = Anorm; % upper bound on the norm of the measurement operator A*G
        param_HSI.reweight_alpha_ff = 0.9;
        param_HSI.reweight_abs_of_max = 0.005;
        param_HSI.use_reweight_steps = 1;
        param_HSI.total_reweights = 30;
        param_HSI.use_reweight_eps = 0;
        
        disp('Split L21 + Nuclear + wavelets')
        disp('-----------------------------------------')
        
        % spectral tesselation (non-overlapping)
        rg_c = domain_decomposition(Qc2, channels(end));
        cell_c_chunks = cell(Qc2, 1);
        y_spmd = cell(Qc2, 1);
        epsilon_spmd = cell(Qc2, 1);
        aW_spmd = cell(Qc2, 1);
        W_spmd = cell(Qc2, 1);
        G_spmd = cell(Qc2, 1);
        
        for i = 1:Qc2
            cell_c_chunks{i} = rg_c(i, 1):rg_c(i, 2);
            y_spmd{i} = y(cell_c_chunks{i});
            epsilon_spmd{i} = epsilons(cell_c_chunks{i});
            aW_spmd{i} = aW(cell_c_chunks{i});
            W_spmd{i} = W(cell_c_chunks{i});
            G_spmd{i} = G(cell_c_chunks{i});
        end
        
        param_HSI.ind = 0;
        
        clear y epsilon aW W G
        
        switch parallel_version            
            case 'standard' 
                % reference version, overlap between the facets underlying the nuclear norms
                [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                    facetHyperSARA(y_spmd, epsilon_spmd, ...
                    A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, Qc2, ...
                    wlt_basis, L, nlevel, cell_c_chunks, channels(end), 'whatever.mat'); % [10/10/2019] ok
                
            case 'weighted' 
                % (piece-wise constant weights, following Audrey's suggestion)
                [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                    facetHyperSARA_weighted(y_spmd, epsilon_spmd, ...
                    A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, Qc2, ...
                    wlt_basis, L, nlevel, cell_c_chunks, channels(end), 'whatever.mat'); % [10/10/2019] ok
                
            case 'cst' 
                % cst overlap for the facets undelrying the nuclear norms (d < overlap from sdwt2)
                % d = (power(2, nlevel)-1)*(max(L(:))-2)/2; % use of db8
                [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                    facetHyperSARA_cst_overlap(y_spmd, epsilon_spmd, ...
                    A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, Qc2, ...
                    wlt_basis, L, nlevel, cell_c_chunks, channels(end), d, 'whatever.mat'); % [10/10/2019] ok
                
            case 'cst_weighted' 
                % same as spmd_sct, weight correction (apodization window in this case)
                [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                    facetHyperSARA_cst_overlap_weighted(y_spmd, epsilon_spmd, ...
                    A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, Qc2, ...
                    wlt_basis, L, nlevel, cell_c_chunks, channels(end), d, window_type, 'whatever.mat') % [10/10/2019] ok
                
                %[xsol,param_HSI,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter] = ...
                %    facetHyperSARA_cst_overlap_weighted_real_data(y_spmd, epsilon_spmd, ...
                %    A, At, aW_spmd, G_spmd, W_spmd, param_HSI, Qx, Qy, Qc2, wlt_basis, L, ...
                %    nlevel, cell_c_chunks, channels(end), d, window_type, 'whatever.mat'); % [10/10/2019] ok basic debugging
                
            case 'standard2' 
                % alternative implementation (gather primal variables and data on the same nodes)
                % gather image and data on the same nodes (extra communications compared to spmd4 for reweigthing and monitoring variables)
                [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                    facetHyperSARA2(y_spmd, epsilon_spmd, ...
                    A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, Qc2, ...
                    wlt_basis, L, nlevel, cell_c_chunks, channels(end), 'whatever.mat'); % [10/10/2019] ok
            
            case 'no_overlap' 
                % same as spmd4, but no overlap for the facets on which the nuclear norms are taken
                [xsol,param,epsilon,t,rel_fval,nuclear,l21,norm_res_out,end_iter,snr_x,snr_x_average] = ...
                    facetHyperSARA_no_overlap(y_spmd, epsilon_spmd, ...
                    A, At, aW_spmd, G_spmd, W_spmd, param_HSI, X0, Qx, Qy, Qc2, ...
                    wlt_basis, L, nlevel, cell_c_chunks, channels(end), 'whatever.mat'); % [10/10/2019] ok
            
            otherwise
                error('Unknown solver version.')
        end
        
        c = size(xsol,3);
        sol = reshape(xsol(:),numel(xsol(:))/c,c);
        snr_x = 20*log10(norm(X0(:))/norm(X0(:)-sol(:)))
        psnrh = zeros(c,1);
        for i = channels
            psnrh(i) = 20*log10(norm(X0(:,i))/norm(X0(:,i)-sol(:,i)));
        end
        snr_x_average = mean(psnrh)
        time_iter_average = mean(end_iter)
        
        mkdir('results/')
        save(fullfile(results_path, results_name),'-v7.3','xsol', 'sol', 'X0', 'SNR', 'SNR_average', 'res');
        fitswrite(xsol,strcat(results_path, 'x_hyperSARA_', alg_version, '_', parallel_version, '_Qx=', num2str(Qx), '_Qy=', num2str(Qy), '_Qc=', num2str(Qc), '.fits'))
        
    end
    
end
