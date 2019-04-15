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
if gaussian_UV
    
    % samplig pattern parameters
    % options 'gaussian', 'file', 'gaussian+large-holes', 'file+undersample', ''gaussian+missing-pixels'
    % sampling_pattern = 'file+undersample';
    % sampling_pattern = 'gaussian';
    sampling_pattern = 'gaussian+large-holes';
    %  sampling_pattern = 'gaussian+missing-pixels';
    % sampling_pattern = 'file';
    
    % sparam.file_name = '/Volumes/Data/MeasSets/meerkat2h.ar.uvw.dat'; % file name for uv coverage
    sparam.file_name = '/Volumes/Data/MeasSets/ska.2h.ar.uvw.dat'; % file name for uv coverage
    sparam.p = percentage(k); % number of measurements as proportion of number of pixels to recover
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
    
    % generate the sampling pattern
    sparam.N = N; % number of pixels in the image
    sparam.Nox = ox*Nx; % number of pixels in the image
    sparam.Noy = oy*Ny; % number of pixels in the image
    
    if generate_uv
        [u, v] = util_gen_sampling_pattern(sampling_pattern, sparam);
        if save_uv
            save(['hypersara-sdwt2/simulated_data/data/uv_60_HI_final_' num2str(percentage(k)) '.mat'],'-v7.3', 'u', 'v');
        end
    elseif load_uv
        load(['hypersara-sdwt2/simulated_data/data/uv_60_HI_final_' num2str(percentage(k)) '.mat']);
        disp('coverage loaded successfully ;)')
    end
    
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
    
    if ~simple_test_1B
        u1 = u1/2;
        v1 = v1/2;
    end
    
    %scatter(u1,v1,'r.')
end

if realistic_UV
    param_data.cov_type = 'random'; 'vlaa';
    param_data.T = 100; % 200
    param_data.na = 25;
    
    uv(1,:) = rand(1,2) ;                                   %
    for alpha = 2:param_data.na                             %
        uv_ = rand(1,2) ;                                       %
        while ismember(uv_,uv,'rows')                           %
            uv_ = rand(1,2) ;                                       %
        end                                                     %
        uv(alpha,:) = uv_ ;                                     %
    end                                                     %
    Pos_ant = 1e06*[uv, zeros(param_data.na,1)] ;           %
    save_cov_file = ['gen_uv_tracks/rand_pos.mat'];         %
    
    %load('ant_vla_pos.mat');
    %load('rand_pos_na50_nstep10_seed1.mat');
    param_data.na = max(size(Pos_ant)) ;
    param_data.M = param_data.T*param_data.na*(param_data.na-1)/2 ;
    param_data.ant_pair = param_data.na*(param_data.na-1)/2 ;
    
    % -------------------------------------------------------------------------
    h = linspace(-10, 10, param_data.T)*pi/12;% hour angle range of +/-  hours -4, 4
    % % % h = linspace(-3, 3, param_data.T)*pi/12;% hour angle range of +/-  hours
    dec = (pi/180)*(40);  % Cas A is 56.4
    lat = (pi/180)*(38. + 59./60. + 48./3600.);  % College Park
    % position de reference x00 dans le plan x-y
    x00 = [mean(Pos_ant(:,1)), mean(Pos_ant(:,2))] ;
    % -------------------------------------------------------------------------
    
    if generate_uv
        [U,V,~] = generate_uv_cov_antennas(Pos_ant,x00,h,lat,dec, param_data.T) ;
        if save_uv
            save(['hypersara-sdwt2/simulated_data/data/uv_60_HI_final_' num2str(percentage(k)) '.mat'],'-v7.3', 'U', 'V');
        end
    elseif load_uv
        load(['hypersara-sdwt2/simulated_data/data/uv_60_HI_final_' num2str(percentage(k)) '.mat']);
        disp('coverage loaded successfully ;)')
    end
    
    
    uab = cell(1,param_data.T) ;
    vab = cell(1,param_data.T) ;
    for t=1:param_data.T
        uab{t} = zeros(param_data.ant_pair, 1) ;
        vab{t} = zeros(param_data.ant_pair, 1) ;
        i = 1 ;
        for alp = 1:param_data.na
            for bet = alp+1:param_data.na
                uab{t}(i) = U{alp}(t)- U{bet}(t) ;
                vab{t}(i) = V{alp}(t)- V{bet}(t) ;
                i=i+1 ;
            end
        end
    end
    
    uu = cell2mat(uab) ; vv = cell2mat(vab) ;
    maxuv = max(max(abs(uu(:))), max(abs(vv(:)))) ;
    uuu = uu/maxuv*pi*0.95 ; vvv = vv/maxuv*pi*0.95 ;
    
    %
    % u = cell2mat(U) ; u = u(:) ;
    % v = cell2mat(V) ; v = v(:) ;
    % maxuv = max(max(abs(u)), max(abs(v))) ;
    % u = u/maxuv*pi ;
    % v = v/maxuv*pi ;
    
    u1 = uuu(:);
    v1 = vvv(:) ;
    
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
    
    %figure, scatter(u1,v1,'r.');
    save('hypersara-sdwt2/results/uv.mat','u1','v1'); 
 
    if ~generate_simple_image
        u1 = u1/2;
        v1 = v1/2;
    end
    
end

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
    tstart = tic;
    
    [A, At, G{i}, W{i}, Gw{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
end
%figure(5), scatter(uw{i},vw{i},'r.')
%% Save data
if save_data
    save('hypersara-sdwt2/simulated_data/data/data.mat','-v7.3', 'G', 'W', 'aW');
end

if save_full_operator && exist('Gw','var')
    save('hypersara-sdwt2/simulated_data/data/Gw.mat','-v7.3', 'Gw');
end

%% Free memory
if free_memory
    clear u v u1 v1 uw vw aWw nW nWw;
end
