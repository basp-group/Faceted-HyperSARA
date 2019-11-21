addpath ../lib/utils/
addpath ../fouRed/
addpath ../lib/operators
addpath ../lib/nufft

%% real data generation
extract_raw_data = false;
ind = 6;
save_data = true;
%visibility_file_name = '/home/h019/aa16/S-band/CYG-S.mat';
visibility_file_name = {'5','7'};
configuration_file_name = {'C','A'};
param_real_data.image_size_Nx = 2560;
param_real_data.image_size_Ny = 1536;
%param_real_data.pixel_size = 0.3; % in arcsec

%% config parameters
Nx = param_real_data.image_size_Nx;
Ny = param_real_data.image_size_Ny;
N = Nx * Ny;

ox = 2; % oversampling factors for nufft
oy = 2; % oversampling factors for nufft
Kx = 7; % number of neighbours for nufft
Ky = 7; % number of neighbours for nufft

%% preconditioning parameters
param_precond.N = N; % number of pixels in the image
param_precond.Nox = ox*Nx; % number of pixels in the image
param_precond.Noy = oy*Ny; % number of pixels in the image
param_precond.gen_uniform_weight_matrix = 1; %set weighting type
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

%% load real_data;
if extract_raw_data
    [y, uw1, vw1, nWw, f, time_, pos] = util_load_real_vis_data(visibility_file_name,configuration_file_name, ind);
    [uw, vw, pixel_size] = scale_uv(visibility_file_name,configuration_file_name, uw1, vw1);
    
    save(['CYG_data_raw_ind=' num2str(ind) '.mat'],'-v7.3', 'y', 'uw', 'vw', 'nWw', 'f', 'time_', 'pos', 'pixel_size');
else
    load(['/lustre/home/shared/sc004/CYG_data_raw_ind=' num2str(ind) '.mat']);
end

figure(1),
for i = length(uw)
    scatter(uw{i},vw{i},'.');hold on;
end


%%
% ch = [1 : size(y,2)];
%ch = [1 : 16];
ch = [1:32];
ii = 1;

for l = 1:length(ch)
    
    i = ch(l)
    
    %% compute weights
    [aWw] = util_gen_preconditioning_matrix(uw{i}, vw{i}, param_precond);
    
    % set the blocks structure
    if param_block_structure.use_manual_partitioning == 1
        param_block.size = 90000; length(uw{i});
        param_block.snapshot = 0;
        param_block.pos = pos{i};
        out_block = util_time_based_block_sp_ar(uw{i},time_{i},param_block);
        partition = out_block.partition;
        partition1 = [sum(partition(1:2)) sum(partition(3:end))]
        param_block_structure.partition = partition1;
    end
    
    % set the blocks structure
    [u, v, ~, uvidx, aW{ii}, nW] = util_gen_block_structure(uw{i}, vw{i}, aWw, nWw{i}, param_block_structure);
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    tstart = tic;
    
    % compute A & At with Kx = Ky = 7 
    [A, At, ~, ~] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW); % Kx = Ky = 7
    
    % load the G matrices and the estimated vis.
    if i < 17
    dataStruct=     load(['/lustre/home/shared/sc004/monoCalibratedWideband4Abdullah/WB5-' num2str(i) '/G-MODELVIS-CELL.mat']);
    else
     dataStruct=load(['/lustre/home/shared/sc004/monoCalibratedWideband4Abdullah/WB7-' num2str(i) '/G-MODELVIS-CELL.mat']);
    end
    Gw = [dataStruct.G{2}; dataStruct.G{1}]; %concatenate the configurations C first and then A
    y_est_full = [dataStruct.VisEstSliced{2}; dataStruct.VisEstSliced{1}]; %concatenate the configurations C first and then A
    clear dataStruct;
    disp('Gw')
    size(Gw)

    % compute the W matrices for all the blocks with the correct convolunatinal kernel size
    %Ky = sqrt(nnz(Gw(1,:))); % this size will apply for all the blocks, is this correct Arwa?
    %Kx = sqrt(nnz(Gw(1,:)));
%%%%%%%% ARWA : %    [~, ~, ~, W{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);

    G{ii} = cell(length(u),1);
    yb{ii} = cell(length(u),1);
    y_est{ii} = cell(length(u),1);
    W{ii}  =cell(length(u),1); %%%%%%%% ARWA
    for j = 1 : length(u)
        G{ii}{j} = Gw(uvidx{j},:);
        W{ii}{j}= G{ii}{j}' * ones(size(G{ii}{j}, 1), 1) ~= 0;

        G{ii}{j} = G{ii}{j}(:, W{ii}{j});
        y_est{ii}{j} = y_est_full(uvidx{j}); 
        yb{ii}{j} = y{i}(uvidx{j});
%         eps_b{ii}{j} = norm(yb{ii}{j} - y_est{ii}{j});
        res{ii}{j} = yb{ii}{j} - y_est{ii}{j};
    end
    clear Gw y_est_full y_est;

%     ii = ii + 1;
    
    save(['/lustre/home/shared/sc004/cyg_data_sub/res_2b_ch', num2str(i),'_ind=' num2str(ind) '.mat'],'-v7.3', 'res');
    save(['/lustre/home/shared/sc004/cyg_data_sub/CYG_data_cal_2b_ch', num2str(i),'_ind=', num2str(ind),'.mat'],'-v7.3', 'G', 'W', 'aW', 'yb');
end
% save(['res_ind=' num2str(ind) '.mat'],'-v7.3', 'res');

%% Save data
% if save_data
%     save('CYG_data_cal_2b_ch32_ind=6.mat','-v7.3', 'G', 'aW', 'yb');
%     save('CYG_y.mat','-v7.3', 'y');
% end

%%
% if save_full_operator && exist('Gw','var')
%     save('CYG_Gw.mat','-v7.3', 'Gw');
% end
% 
% 
% %% Free memory
% if free_memory
%     clear y u v uv_mat uv_mat1 uvidx uw vw ant1 ant2 aWw nW nWw out_block;
% end
