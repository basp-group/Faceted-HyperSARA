%% config parameters
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

%% samplig pattern parameters
% options 'gaussian', 'file', 'gaussian+large-holes', 'file+undersample', ''gaussian+missing-pixels'
% sampling_pattern = 'file+undersample';
% sampling_pattern = 'gaussian';
sampling_pattern = 'gaussian+large-holes';
%  sampling_pattern = 'gaussian+missing-pixels';
% sampling_pattern = 'file';

% sparam.file_name = '/Volumes/Data/MeasSets/meerkat2h.ar.uvw.dat'; % file name for uv coverage
sparam.file_name = '/Volumes/Data/MeasSets/ska.2h.ar.uvw.dat'; % file name for uv coverage
sparam.p = percentage; % number of measurements as proportion of number of pixels to recover
sparam.hole_number = 8000; % number of holes to introduce for 'gaussian+large-holes'
sparam.hole_prob = 0.05; % probability of single pixel hole for 'gaussian+missing-pixels'
sparam.hole_size = pi/60; % size of the missing frequency data
sparam.fpartition = [pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [0, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.25*pi, 0, 0.25*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.3*pi, -0.15*pi, 0, 0.15*pi, 0.3*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
% sparam.fpartition = [-0.35*pi, -0.25*pi, -0.15*pi, 0, 0.15*pi, 0.25*pi, 0.35*pi, pi]; % partition (symetrically) of the data to nodes (frequency ranges)
sparam.sigma = pi/3; % variance of the gaussion over continous frequency
sparam.sigma_holes = pi/3; % variance of the gaussion for the holes

%% generate the sampling pattern
sparam.N = N; % number of pixels in the image
sparam.Nox = ox*Nx; % number of pixels in the image
sparam.Noy = oy*Ny; % number of pixels in the image


%%
if generate_uv
    if ~ usingCalibrationKernels
        [u, v] = util_gen_sampling_pattern(sampling_pattern, sparam);
    else
        % to be modified, load realistic uv-coverage
        dl = 1.1;
%         [u, v, na] = generate_uv_coverage(T, hrs, dl, 'vlaa');
        load('data/test_DR.mat', 'uvw');
        u = uvw(:,1);
        v = uvw(:,2);
        u = {u(:)};
        v = {v(:)};
        na = 27;
        M = na*(na-1)/2;
    end
    
    if save_data
        save(['./simulated_data/data/uv_' num2str(c) '_HI_final_' num2str(percentage(k)) '.mat'],'-v7.3', 'u', 'v');
    end
else
    load(['./simulated_data/data/uv_' num2str(c) '_HI_final_' num2str(percentage(k)) '.mat']);
    disp('coverage loaded successfully ;)')
end

% % remove all the visibilities outside the range [-pi, pi]
% r = sqrt(u{1}.^2 + v{1}.^2);
% size(r(r>pi))
% bmax = max(r)
% 
% u1 = u{1};
% v1 = v{1};
% [mm] = find(r > pi);
% u1(mm) = 0;
% v1(mm) = 0;
% 
% u1 = u1(u1 ~= 0);
% v1 = v1(v1 ~= 0);
% 
% r = sqrt(u1.^2 + v1.^2);
% size(r(r>pi))
% bmax = max(r)
% 
% u1 = u1/2;
% v1 = v1/2;

% rescale visibilities
bmax = max(sqrt(u{1}.^2 + v{1}.^2));

% setting imaging params
v1 = v{1}*pi/(bmax*dl); % dl = 1.1
u1 = u{1}*pi/(bmax*dl);

%scatter(u1,v1,'r.')

%%

if usingCalibrationKernels
    S = 3;
    S2 = S^2;
    cs = floor(S2/2) + 1; % center of the spatial support
    D = 1e-3*(randn(S2, na, T) + 1i*randn(S2, na, T));
    D(cs, :) = 1;
end

for i = 1:length(ch)
    
    uw{i} = (f(i)/f(1)) * u1;
    vw{i} = (f(i)/f(1)) * v1;
 
%     % use the absolute values to speed up the search
%     Gw_a = abs(Gw{i});
% 
%     b_l = length(uw{i});
%     % check if eack line is entirely zero
%     W{i} = Gw_a' * ones(b_l, 1) ~= 0;
% 
%     % store only what we need from G
%     G{i} = Gw{i}(:, W{i});
    
    %% compute uniform weights (sampling density) for the preconditioning
    [aWw] = util_gen_preconditioning_matrix(uw{i}, vw{i}, param_precond);
    
    % set the weighting matrix, these are the natural weights for real data
    nWw = ones(length(uw{i}), 1);
    
    % set the blocks structure
    [u, v, ~, uvidx, aW{i}, nW] = util_gen_block_structure(uw{i}, vw{i}, aWw, nWw, param_block_structure);
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    %tstart = tic;
    if usingBlocking
        if usingCalibrationKernels
            [A, At, G{i}, W{i}] = op_p_nufft_calib([v u], D, [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], S, nW);
        else
            [A, At, G{i}, W{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
        end
    else
        if usingCalibrationKernels
            [A, At, G{i}, ~] = op_nufft_calib([vw{i} uw{i}], D, [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], S);
        else
            [A, At, G{i}, ~] = op_nufft([vw{i} uw{i}], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);
        end 
    end
end

%% Free memory
clear u v u1 v1 uw vw aWw nW nWw r antennas na mm bmax uvidx Gw_a W b_l;
