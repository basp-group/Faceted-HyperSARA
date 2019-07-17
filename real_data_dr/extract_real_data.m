%% real data generation
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
Kx = 8; % number of neighbours for nufft
Ky = 8; % number of neighbours for nufft

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
if extract_raw_data % MODIFY READING HERE (MATFILE)
    [y, uw1, vw1, nWw, f, time_, pos] = util_load_real_vis_data(visibility_file_name,configuration_file_name, ind);
    [uw, vw, pixel_size] = scale_uv(visibility_file_name,configuration_file_name, uw1, vw1);
    
    save(['CYG_data_raw_ind=' num2str(ind) '.mat'],'-v7.3', 'y', 'uw', 'vw', 'nWw', 'f', 'time_', 'pos', 'pixel_size');
else
    load(['CYG_data_raw_ind=' num2str(ind) '.mat']);
end

figure(1),
for i = length(uw)
    scatter(uw{i},vw{i},'.');hold on;
end


%%
ch = [1 : size(y,2)];
for i = ch
    
    i
    
    %% compute weights
    [aWw] = util_gen_preconditioning_matrix(uw{i}, vw{i}, param_precond);
    
    % set the blocks structure
    if param_block_structure.use_manual_partitioning == 1
        param_block.size = 90000; length(uw{i});
        param_block.snapshot = 0;
        param_block.pos = pos{i};
        out_block = util_time_based_block_sp_ar(uw{i},time_{i},param_block);
        partition = out_block.partition
        param_block_structure.partition = partition;
    end
    
    % set the blocks structure
    [u, v, ~, uvidx, aW{i}, nW] = util_gen_block_structure(uw{i}, vw{i}, aWw, nWw{i}, param_block_structure); % CHANGE DEF. OF OPERATOR
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    tstart = tic;
   
    %     for k = 1 : length(nW)
    %         nW{k} = ones(size(nW{k}));
    %     end
    
    %[A, At, G{i}, W{i}, Gw{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
    [A, At, G{i}, W{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
    
    yb{i} = cell(length(u),1);
    for j = 1 : length(u) 
        yb{i}{j} = y{i}(uvidx{j});
    end
    
end

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
if free_memory
    clear y u v uv_mat uv_mat1 uvidx uw vw ant1 ant2 aWw nW nWw out_block;
end
