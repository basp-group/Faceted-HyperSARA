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
% param_data.cov_type = 'vlaa';
% hrs = 6;

% [u1, v1, na, antennas] = generate_uv_coverage2(T, hrs, param_data.cov_type);
[u, v] = util_gen_sampling_pattern(sampling_pattern, sparam);

% remove all the visibilities outside the range [-pi, pi]
r = sqrt(u{1}.^2 + v{1}.^2);
size(r(r>pi))
bmax = max(r)

u1 = u{1};
v1 = v{1};
[mm] = find(r > pi);
u1(mm) = 0;
v1(mm) = 0;

u1 = u1(u1 ~= 0);
v1 = v1(v1 ~= 0);

r = sqrt(u1.^2 + v1.^2);
size(r(r>pi))
bmax = max(r)

u1 = u1/2;
v1 = v1/2;

%figure(11), scatter(u1,v1,'r.'); %hold on; scatter(-u1,v1,'b.');

%remove all the visibilities outside the range [-pi, pi]
% r = sqrt(u1.^2 + v1.^2);
% size(r(r>pi))
% bmax = max(r)
% 
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
% %figure(2), scatter(u1,-v1,'b.');
% 
% % mkdir('results')
% % save('./results/uv.mat','u1','v1');
% 
% if ~generate_simple_image
%     u1 = u1/2;
%     v1 = v1/2;
% end
% 
% size(u1)

%%
for i = 1:nChannels
    
    uw{i} = (f(i)/f(1)) * u1;
    vw{i} = (f(i)/f(1)) * v1;
    
    %% compute uniform weights (sampling density) for the preconditioning
    [aWw] = util_gen_preconditioning_matrix(uw{i}, vw{i}, param_precond);
    
    % set the weighting matrix, these are the natural weights for real data
    nWw = ones(length(uw{i}), 1);
    
    % set the blocks structure
    [u, v, ~, uvidx, aW{i}, nW] = util_gen_block_structure(uw{i}, vw{i}, aWw, nWw, param_block_structure);
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    %tstart = tic;
   
    [A, At, G{i}, W{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW); 
    %[A, At, G{i}, W{i}, Gw{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
end
%figure, scatter(uw{i},vw{i},'r.');

%% Save data
if save_data
    mkdir('/simulated_data/data')
    save('/simulated_data/data/data.mat','-v7.3', 'G', 'W', 'aW');
end

if save_full_operator && exist('Gw','var')
    mkdir('/simulated_data/data')
    save('/simulated_data/data/Gw.mat','-v7.3', 'Gw');
end

%% Free memory
if free_memory
    clear u v u1 v1 uw vw aWw nW nWw r antennas na mm bmax uvidx f;
end
