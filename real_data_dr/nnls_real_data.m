function nnls_real_data(chInd)
fprintf('Channel number %d', chInd);

%% Real data extraction
addpath ../lib/utils/
addpath ../fouRed/
addpath ../lib/operators
addpath ../lib/nufft
% addpath ../data_mnras_dr

param_real_data.image_size_Nx = 2560;
param_real_data.image_size_Ny = 1536;
% concatenation) for the real dataset considered
nBlocks = 2;        % number of data blocks (needs to be known beforehand,
% quite restrictive here), change l.70 accordingly

%% Config parameters
Nx = param_real_data.image_size_Nx;
Ny = param_real_data.image_size_Ny;
N = Nx * Ny;
ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

new_file_y = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_y.mat');
new_file_u = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_u.mat');
new_file_v = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_v.mat');
new_file_nW = matfile('/lustre/home/shared/sc004/cyg_data_2b_dr/CYG_2b_nW.mat');

%% Estimate epsilon with NNLS on each data block
    
    % set up res_nnls to be saved on disk
res_nnls = cell(1, 1);
for l = 1:numel(res_nnls)
    res_nnls{l} = cell(nBlocks, 1);
end
save(['CYG_2b_res_nnls=', num2str(chInd),'.mat'], 'res_nnls', '-v7.3');
clear res_nnls

% parameter NNLS
param_nnls.verbose = 2;       % print log or not
param_nnls.rel_obj = 1e-10;    % stopping criterion
param_nnls.max_iter = 1000;     % max number of iterations 1000
param_nnls.sol_steps = [inf]; % saves images at the given iterations
param_nnls.beta = 1;

new_file_res = matfile(['CYG_2b_res_nnls=', num2str(chInd),'.mat'], 'Writable', true);

% define operators
parpool(2)
[A, At, ~, ~] = op_nufft([0, 0], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);

% solve NNLS per block / estimate epsilon / reduce data
y_tmp = new_file_y.y(chInd,1);
u_tmp = new_file_u.u(chInd,1);
v_tmp = new_file_v.v(chInd,1);
nW_tmp = new_file_nW.nW(chInd,1);
res_tmp = new_file_res.res_nnls(1,1);
y_tmp = y_tmp{1};
u_tmp = u_tmp{1};
v_tmp = v_tmp{1};
nW_tmp = nW_tmp{1};
res_tmp = res_tmp{1};

parfor j = 1:nBlocks
    [~, ~, G, W1] = op_p_nufft([v_tmp(j) u_tmp(j)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW_tmp(j));
    [sol, ~] = fb_nnls_blocks(y_tmp{j}, A, At, G{1}(:, W1{1}), W1{1}, param_nnls); % check format of "sol" (matrix)
    res_tmp{j} = y_tmp{j} - G{1}*A(sol); % check computation of the residual
end
new_file_res.res_nnls(1,1) = {res_tmp};

fprintf('NNLS of channel %d is finished \n', chInd);
end
