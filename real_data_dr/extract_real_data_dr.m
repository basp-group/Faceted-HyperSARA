%% real data generation
addpath ../Dimensionality-reduced-hyper-SARA-sdwt/lib/utils/
addpath ../Dimensionality-reduced-hyper-SARA-sdwt/fouRed/

%visibility_file_name = '/home/h019/aa16/S-band/CYG-S.mat';
visibility_file_name = 'CYG_data_raw_ind=';
param_real_data.image_size_Nx = 2560;
param_real_data.image_size_Ny = 1536;
nSpw = 2;
nChannels = 2*nSpw*nSpw; % total number of channels for the real dataset 
%                          considered
klargestpercent = 20;
extract_raw_data = true;
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

%% Fourier reduction parameters
param_fouRed.enable_klargestpercent = 1;
param_fouRed.klargestpercent = klargestpercent;
param_fouRed.enable_estimatethreshold = 0;
param_fouRed.gamma = 3;             % By using threshold estimation, the optimal theshold reads as gamma * sigma / ||x||_2
param_fouRed.diagthresholdepsilon = 1e-10; 
param_fouRed.covmatfileexists = 0;
param_fouRed.covmatfile = 'covariancemat.mat';
param_fouRed.fastCov = 1;

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
%     [y, uw1, vw1, nWw, f, time_, pos] = util_load_real_vis_data(visibility_file_name,configuration_file_name, ind);
%     [uw, vw, pixel_size] = scale_uv(visibility_file_name,configuration_file_name, uw1, vw1);
%     save(['CYG_data_raw_ind=' num2str(ind) '.mat'],'-v7.3', 'y', 'uw', 'vw', 'nWw', 'f', 'time_', 'pos', 'pixel_size');

%     % solution 1: 32 independent mat files
%     [y, uw, vw, nWw, f, time_, pos] = util_load_real_vis_data_dr(visibility_file_name, nSpw, ind);
%     save(['CYG_data_raw_dr_ind=' num2str(ind) '.mat'], '-v7.3', 'y', 'uw', 'vw', 'nWw', 'f', 'time_', 'pos');
    
    % solution 2: 1 dataset for each variable of interest (see for the 
    % blocking strategy to be applied here)
%     y = cell(32, 1);
%     for l = 1:2*nSpw
%         y{l} = [];
%         for n = 1:nSpw
%             filename = [visibility_file_name, num2str(n), '.mat'];
%             file = matfile(filename);  
%             y{l} = [y{l}; cell2mat(file.y(1,l))];
%         end
%     end
    % or 
%     for l = 1:2*nSpw
%         y{l} = cell(16, 1);
%         for n = 1:nSpw
%             filename = [visibility_file_name, num2str(n), '.mat'];
%             file = matfile(filename);  
%             y{l}{n} = cell2mat(file.y(1,l));
%         end
%     end
    % or (probably the best option so far) 'y', 'uw', 'vw', 'nWw', 'f', 'time_', 'pos'
    k = 1;
    y = cell(2*nSpw*nSpw, 1);
    for l = 1:2*nSpw
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            y{k} = cell2mat(file.y(1,l));
            k = k+1;
        end
    end
    save('CYG_y.mat', 'y', '-v7.3');
    clear y;
    
    k = 1;
    nWw = cell(2*nSpw*nSpw, 1);
    for l = 1:2*nSpw
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            nWw{k} = cell2mat(file.nWw(1,l));
            k = k+1;
        end
    end
    save('CYG_nWw.mat', 'nWw', '-v7.3');
    clear y;
    
    k = 1;
    uw = cell(2*nSpw*nSpw, 1);
    for l = 1:2*nSpw
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            uw{k} = cell2mat(file.uw(1,l));
            k = k+1;
        end
    end
    save('CYG_uw.mat', 'uw', '-v7.3');
    clear uw
    
    k = 1;
    vw = cell(2*nSpw*nSpw, 1);
    for l = 1:2*nSpw
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            vw{k} = cell2mat(file.vw(1,l));
            k = k+1;
        end
    end
    save('CYG_vw.mat', 'vw', '-v7.3');
    clear vw
    
    k = 1;
    f = cell(2*nSpw*nSpw, 1);
    for l = 1:2*nSpw
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            f{k} = cell2mat(file.f(1,l));
            k = k+1;
        end
    end
    save('CYG_f.mat', 'f', '-v7.3');
    clear f
    
    k = 1;
    time_ = cell(2*nSpw*nSpw, 1);
    for l = 1:2*nSpw
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            time_{k} = cell2mat(file.time_(1,l)); % directly concatenated, or not?
            k = k+1;
        end
    end
    save('CYG_time.mat', 'time_', '-v7.3');
    clear time_
    
    k = 1;
    pos = cell(2*nSpw*nSpw, 1);
    for l = 1:2*nSpw
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            pos{k} = cell2mat(file.pos(1,l));
            k = k+1;
        end
    end
    save('CYG_pos.mat', 'pos', '-v7.3');
    clear pos
    % see what needs to be finalized in terms of dimensionality reduction
    % idem for the reamining variables: H, ...
else
    %load(['CYG_data_raw_ind=' num2str(ind) '.mat']);
    load(['CYG_data_raw_dr_ind=' num2str(ind) '.mat']); % load only the variables needed at one point (oterhwise, quite large!)
end

% figure(1),
% for i = length(uw)
%     scatter(uw{i},vw{i},'.');hold on;
% end


%% adapt to DR (see codes from Ming / ask Ming directly if necessary)
%%
ch = [1 : numel(y)];
% output from the definition of the Fourier reduction scheme
H = cell(nChannels, 1);  % holographic matrices
yT = cell(nChannels, 1); % reduced data
T = cell(nChannels, 1);  % preconditioning matrix (inverse of the singular values)
Wm = cell(nChannels, 1); % masks
aW = cell(nChannels, 1); % 1./T (shall I keep this variable?)

u = cell(nChannels, 1);
v = cell(nChannels, 1);
uvidx = cell(nChannels, 1);

for i = ch
    
    i
    
    file_y = matfile('CYG_y.mat');
    file_pos = matfile('CYG_pos.mat');
    file_uw = matfile('CYG_uw.mat'); 
    file_vw = matfile('CYG_vw.mat');
    file_time_ = matfile('CYG_time.mat');
    file_nWw = matfile('CYG_nWw.mat');
    
    %% compute weights
    uw = cell2mat(file_uw.uw(i, 1));
    vw = cell2mat(file_vw.vw(i, 1));
    pos = cell2mat(file_pos.pos(i, 1));
    time_ = cell2mat(file_time.time_(i, 1));
    %[aWw] = util_gen_preconditioning_matrix(uw, vw, param_precond); % see if still necessary
    
    % set the blocks structure
    if param_block_structure.use_manual_partitioning == 1
        param_block.size = 90000; length(uw);
        param_block.snapshot = 0;
        param_block.pos = pos;
        out_block = util_time_based_block_sp_ar(uw,time_,param_block);
        param_block_structure.partition = out_block.partition;
    end
    
    u{i} = [];
    v{i} = [];
    nW{i} = [];
    uvidx{i} = [];
    for k = 1:numel(uw)
        % set the blocks structure (per frequency) [is aW{i} necessary?]
        [u1, v1, ~, uvidx1, aW, nW1] = util_gen_block_structure(uw{k}, vw{k}, aWw{k}, nWw{k}, param_block_structure);
        
        % reorder the blocks along a single dimension (remove intermediate level of the frequencies)
        u{i} = [u{i}; u1];
        v{i} = [v{i}; v1];
        nW{i} = [nW{i}; nW1];
        uvidx{i} = [uvidx{i}; uvidx1];
    end
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    tstart = tic;
   
    %     for k = 1 : length(nW)
    %         nW{k} = ones(size(nW{k}));
    %     end
    
    %[A, At, G{i}, W{i}, Gw{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
    % do I need to keep all the G matrices? (to be determined later on)
    %[A, At, G{i}, W{i}] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
    [A, At, G, ~] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW);
    
    %% --- Fourier reduction --- %%
    H{i} = cell(numel(G), 1);
    T{i} = cell(numel(G), 1);
    yT{i} = cell(numel(G), 1);
    Wm{i} = cell(numel(G), 1);
    for j = 1:numel(G)
        
        H{i}{j} = (G{j}')*G{j};
        
        % estimate threshold
        dirty2 = norm(operatorPhit(y{i}{j}, G{j}', At) / sqrt(N));
        
        ... define Mask with the appropriate threshold -> W (how is this supposed to be done on real data?)
        % no longer invoke function fourierReduction to reduce lambda functions
        % fast way of matrix probing (using psf)
        dirac2D = zeros(Ny, Nx);
        dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;
        PSF = operatorIpsf(dirac2D, A, At, H{i}{j}, [oy*Ny, ox*Nx]);
        covariancemat = FT2(PSF);
        d_mat = abs(real(covariancemat(:)));
        %
        if param_fouRed.enable_klargestpercent
            Mask = (d_mat >= prctile(d_mat,100-param_fouRed.klargestpercent));
        elseif param_fouRed.enable_estimatethreshold
            % Embed the noise
            noise1 = param_fouRed.sigma_noise * (randn(size(G{j}, 1),1) + 1j * randn(size(Gw{i}, 1), 1));
            rn = FT2(At(G{j}'*noise1));  % Apply F Phi
            th_dirty = param_fouRed.gamma * std(rn(:)) / dirty2;
            fprintf('\nThe estimate threshold using ground truth is %e \n', th);
            Mask = (d_mat >= th_dirty);
        end
        d_mat = d_mat(Mask);
        %
        T{i}{j} = max(param_fouRed.diagthresholdepsilon, d_mat);  % This ensures that inverting the values will not explode in computation
        T{i}{j} = 1./sqrt(T{i}{j});
        Wm{i}{j} = Mask;
        
        % reduce the data (per block)
        yT{i}{j} = dataReduce(y{i}(uvidx{j}), G{j}', At, T{i}{j}, W{i}{j}); % to be adapted
    end
    % missing: definition of epsilon, yT
    %% ---
    
%     % block the data
%     yb{i} = cell(length(u),1);
%     for j = 1 : length(u) 
%         yb{i}{j} = y{i}(uvidx{j});
%     end
    
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
