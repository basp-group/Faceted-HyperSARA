function func_solver_real_data(gamma, chInd)
% chInd = 1;
addpath ../fouRed
addpath ../lib/
addpath ../lib/operators/
addpath ../lib/nufft/
addpath ../lib/utils/
addpath ../lib/CubeHelix/
addpath ../lib/Proximity_operators/code/matlab/indicator/
addpath ../lib/Proximity_operators/code/matlab/multi/
addpath ../sdwt2/
addpath ../src/
addpath ../src/spmd/
addpath ../src/spmd/dr/
addpath ../src/spmd/weighted/

disp(chInd)

compute_Anorm = true;
usingPrecondition = false;
compute_G = false;

flag_algo = 1;

param_real_data.image_size_Nx = 2560; % 2560;
param_real_data.image_size_Ny = 1536; % 1536;
nChannels = 1; % total number of "virtual" channels (i.e., after
% concatenation) for the real dataset considered
nBlocks = 2;        % number of data blocks (needs to be known beforehand,
% quite restrictive here), change l.70 accordingly
% klargestpercent = 20;

%% Config parameters
Nx = param_real_data.image_size_Nx;
Ny = param_real_data.image_size_Ny;
N = Nx * Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

%% Preconditioning parameters
param_precond.N = N;       % number of pixels in the image
param_precond.Nox = ox*Nx; % number of pixels in the image
param_precond.Noy = oy*Ny; % number of pixels in the image
param_precond.gen_uniform_weight_matrix = 1; % weighting type
param_precond.uniform_weight_sub_pixels = 1;

%% block structure
regenerate_block_structure = 1;

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

%% Load data
% % data files in EPFL
% new_file_y = matfile('/Users/ming/workspace/Git/extract_real_data/CYG_data_cal_9b_ch32_ind=6.mat');
% new_file_res = matfile('/Users/ming/workspace/Git/extract_real_data/res_9b_ch32_ind=6.mat');
% data files on cirrus
new_file_y = matfile('/lustre/home/shared/sc004/cyg_data_sub/old/CYG_data_cal_2b_ch32_ind=6.mat');
new_file_res = matfile('/lustre/home/shared/sc004/cyg_data_sub/old/res_2b_ch32_ind=6.mat');
% yb = new_file_y.y(chInd,1);
yb = new_file_y.yb(1,1);
G = new_file_y.G(1,1);
aW = new_file_y.aW(1,1);
res = new_file_res.res(1,1);
for j = 1:nBlocks
    epsilon{1}{j} = norm(res{1}{j});
    W{1}{j} = true(Nx*Ny*ox*oy,1);
end

[A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);

% if compute_G  
%     new_file_u = matfile('/Users/ming/workspace/Git/extract_real_data/CYG_2b_u.mat');
%     new_file_v = matfile('/Users/ming/workspace/Git/extract_real_data/CYG_2b_v.mat');
%     new_file_nW = matfile('/Users/ming/workspace/Git/extract_real_data/CYG_2b_nW.mat');
%     new_file_res = matfile(['/Users/ming/workspace/Git/res_nnls/CYG_2b_res_nnls=', num2str(chInd),'.mat']);
% 
%     % solve NNLS per block / estimate epsilon / reduce data
%     u_tmp = new_file_u.u(chInd,1);
%     v_tmp = new_file_v.v(chInd,1);
%     nW_tmp = new_file_nW.nW(chInd,1);
%     res_tmp = new_file_res.res_nnls(1,1);
%     % y_tmp = y_tmp{1};
%     u_tmp = u_tmp{1};
%     v_tmp = v_tmp{1};
%     nW_tmp = nW_tmp{1};
%     res_tmp = res_tmp{1};
% 
%     for j = 1:nBlocks
%         [aW{1}{j}] = util_gen_preconditioning_matrix(v_tmp{j}, u_tmp{j}, param_precond);
%         epsilon{1}{j,1} = norm(res_tmp{j}, 2);
%     end
%     %     [u, v, ~, uvidx, aW{1}, nW] = util_gen_block_structure(u_tmp, v_tmp, aWw, nW_tmp, param_block_structure);
%     [A, At, G{1}, W{1}] = op_p_nufft([v_tmp u_tmp], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW_tmp);
%     save(['SARA_ch', num2str(chInd),'.mat'], '-v7.3', 'G', 'W', 'aW', 'A','At','epsilon');
% else
%     load(['SARA_ch', num2str(chInd),'.mat'])
% end


%% Compute full measurement operator spectral norm
if compute_Anorm
    F = afclean( @(x) HS_forward_operator_precond_G(x, G, A, aW));
    Ft = afclean( @(y) HS_adjoint_operator_precond_G(y, G, At, aW, Ny, Nx));
    Anorm = pow_method_op(F, Ft, [Ny Nx length(chInd)]);   
    save(['Anorm_2b_per_', num2str(chInd), '.mat'],'-v7.3', 'Anorm');
else
    load(['Anorm_2b_per_', num2str(chInd), '.mat']);
end

Anorm

clear F Ft;

%% Sparsity operator definition
nlevel = 4; % wavelet level
wlt_basis = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'}; % wavelet basis to be used
L = [2*(1:8)'; 0]; % length of the filters (0 corresponding to the 'self' basis)

if flag_algo < 2
    
    [Psi1, Psit1] = op_p_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    P = length(Psi1);

    for k = 1 : P
        f = '@(x_wave) HS_forward_sparsity(x_wave,Psi1{';
        f = sprintf('%s%i},Ny,Nx);', f,k);
        Psi{k} = eval(f);

        b(k) = size(Psit1{k}(zeros(Ny,Nx,1)),1);

        ft = ['@(x) HS_adjoint_sparsity(x,Psit1{' num2str(k) '},b(' num2str(k) '));'];
        Psit{k} = eval(ft);
    end

    %% Full sparsity operator
    [Psiw, Psitw] = op_sp_wlt_basis(wlt_basis, nlevel, Ny, Nx);
    bb = size(Psitw(zeros(Ny,Nx,1)),1);

    Psi_full = @(x_wave) HS_forward_sparsity(x_wave,Psiw,Ny,Nx);
    Psit_full = @(x) HS_adjoint_sparsity(x,Psitw,bb);
    
end

disp('Split L1 + wavelets')
% %% PDFB parameter structure sent to the algorithm
param_pdfb.verbose = 2; % print log or not
param_pdfb.nu1 = 1; % bound on the norm of the operator Psi
param_pdfb.nu2 = Anorm; % bound on the norm of the operator A*G
param_pdfb.gamma = gamma; % convergence parameter L1 (soft th parameter)
param_pdfb.tau = 0.49; % forward descent step size
param_pdfb.rel_obj = 1e-10; % stopping criterion
param_pdfb.max_iter = 10000; % max number of iterations
param_pdfb.lambda0 = 1; % relaxation step for primal update
param_pdfb.lambda1 = 1; % relaxation step for L1 dual update
param_pdfb.lambda2 = 1; % relaxation step for L2 dual update
param_pdfb.omega1 = 1;
param_pdfb.omega2 = 1;
param_pdfb.sol_steps = [inf]; % saves images at the given iterations

param_pdfb.use_proj_elipse_fb = logical(usingPrecondition);
param_pdfb.elipse_proj_max_iter = 20;
param_pdfb.elipse_proj_min_iter = 1;
param_pdfb.elipse_proj_eps = 1e-8; % precision of the projection onto the ellipsoid

param_pdfb.use_adapt_eps = 0; % flag to activate adaptive epsilon (Note that there is no need to use the adaptive strategy on simulations)
param_pdfb.adapt_eps_start = 300; % minimum num of iter before stating adjustment
param_pdfb.adapt_eps_tol_in = 0.99; % tolerance inside the l2 ball
param_pdfb.adapt_eps_tol_out = 1.001; % tolerance outside the l2 ball
param_pdfb.adapt_eps_steps = 100; % min num of iter between consecutive updates
param_pdfb.adapt_eps_rel_obj = 5e-4; % bound on the relative change of the solution
param_pdfb.adapt_eps_change_percentage = 0.5*(sqrt(5)-1); % the weight of the update w.r.t the l2 norm of the residual data

param_pdfb.reweight_alpha = (0.8)^20; % the parameter associated with the weight update equation and decreased after each reweight by percentage defined in the next parameter
param_pdfb.reweight_alpha_ff = 0.8;
param_pdfb.total_reweights = 30; % -1 if you don't want reweighting
param_pdfb.reweight_abs_of_max = Inf; % (reweight_abs_of_max * max) this is assumed true signal and hence will have weights equal to zero => it wont be penalised

param_pdfb.use_reweight_steps = 0; % reweighting by fixed steps
param_pdfb.reweight_step_size = 300; % reweighting step size
param_pdfb.reweight_steps = [5000: param_pdfb.reweight_step_size :10000];
param_pdfb.step_flag = 1;

param_pdfb.use_reweight_eps = 1; % reweighting w.r.t the relative change of the solution
param_pdfb.reweight_max_reweight_itr = param_pdfb.max_iter - param_pdfb.reweight_step_size;
param_pdfb.reweight_rel_obj = 1e-4; % criterion for performing reweighting
param_pdfb.reweight_min_steps_rel_obj = 300; % min num of iter between reweights

param_pdfb.rw_tol = 5000;

% param_pdfb.xstar = fitsread('xsol_sara_ch32_calibrated_5e-6_15rw.fits');
% solvers
mkdir('results/')
[xsol,v1,v2,g,weights1,proj,t_block,reweight_alpha,epsilon,t,rel_fval,l11,norm_res,res] = ...
    pdfb_L11_Adapt_blocks_rw_par_precond_new(yb, epsilon, A, At, aW, G, W, Psi, Psit, param_pdfb, chInd);

save(['results/results_SARA_fouRed_ch', num2str(chInd),'.mat'], '-v7.3', ...
    'xsol', 'epsilon', 't', 'rel_fval', 'l11', 'norm_res', 'res');
    