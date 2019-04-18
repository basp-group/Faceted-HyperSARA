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


%%
param_data.cov_type = 'vlaa';

T = 100;
hrs = 6;

[u1, v1, na, antennas] = generate_uv_coverage2(T, hrs, param_data.cov_type);
u1 = u1(:);
v1 = v1(:);

%figure(11), scatter(u1,v1,'r.'); %hold on; scatter(-u1,v1,'b.');

%remove all the visibilities outside the range [-pi, pi]
r = sqrt(u1.^2 + v1.^2);
size(r(r>pi))
bmax = max(r)

[mm] = find(r > pi);
u1(mm) = 0;
v1(mm) = 0;

u1 = u1(u1 ~= 0);
v1 = v1(v1 ~= 0);

r = sqrt(u1.^2 + v1.^2);
size(r(r>pi))
bmax = max(r)

%figure(2), scatter(u1,-v1,'b.');

% mkdir('results')
% save('./results/uv.mat','u1','v1');

if ~generate_simple_image
    u1 = u1/2;
    v1 = v1/2;
end

size(u1)

%%
for i = ch
    
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
    
    [A, At, G{i}, W{i}, Gw{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
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
    clear u v u1 v1 uw vw aWw nW nWw;
end
