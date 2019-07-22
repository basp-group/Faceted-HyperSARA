%% real data generation
addpath ../lib/utils/
addpath ../fouRed/
addpath ../../real_data_dr

%visibility_file_name = '/home/h019/aa16/S-band/CYG-S.mat';
visibility_file_name = 'CYG_data_raw_ind=';
param_real_data.image_size_Nx = 2560;
param_real_data.image_size_Ny = 1536;
nSpw = 2;
nChannels = 2*nSpw; % total number of channels for the real dataset 
%                     considered
% 2*nSpw*nSpw
klargestpercent = 20;
extract_raw_data = true;
%param_real_data.pixel_size = 0.3; % in arcsec
FT2 = @(x) fftshift(fft2(ifftshift(x)));

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
    % or (probably the best option so far) 'y', 'uw', 'vw', 'nWw', 'f', 'time_', 'pos'
    y = cell(2*nSpw, 1);
    for l = 1:2*nSpw
        y{l} = cell(nSpw, 1);
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            y{l}{n} = cell2mat(file.y(1,l));
        end
    end
    save('CYG_y.mat', 'y', '-v7.3');
    clear y;
    
    nWw = cell(2*nSpw, 1);
    for l = 1:2*nSpw
        nWw{l} = cell(nSpw, 1);
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            nWw{l}{n} = cell2mat(file.nWw(1,l));
        end
    end
    save('CYG_nWw.mat', 'nWw', '-v7.3');
    clear y;
    
    uw = cell(2*nSpw, 1);
    for l = 1:2*nSpw
        uw{l} = cell(nSpw, 1);
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            uw{l}{n} = cell2mat(file.uw(1,l));
        end
    end
    save('CYG_uw.mat', 'uw', '-v7.3');
    clear uw
    
    vw = cell(2*nSpw, 1);
    for l = 1:2*nSpw
        vw{l} = cell(nSpw, 1);
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            vw{l}{l} = cell2mat(file.vw(1,l));
        end
    end
    save('CYG_vw.mat', 'vw', '-v7.3');
    clear vw
    
    f = cell(2*nSpw, 1);
    for l = 1:2*nSpw
        f{l} = cell(nSpw, 1);
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            f{l}{n} = cell2mat(file.f(1,l));
        end
    end
    save('CYG_f.mat', 'f', '-v7.3');
    clear f
    
    time_ = cell(2*nSpw, 1);
    for l = 1:2*nSpw
        time_{l} = cell(nSpw, 1);
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            time_{l}{n} = cell2mat(file.time_(1,l)); % directly concatenated, or not?
        end
    end
    save('CYG_time.mat', 'time_', '-v7.3');
    clear time_
    
    pos = cell(2*nSpw, 1);
    for l = 1:2*nSpw
        pos{l} = cell(nSpw, 1);
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);  
            pos{l}{n} = cell2mat(file.pos(1,l));
        end
    end
    save('CYG_pos.mat', 'pos', '-v7.3');
    clear pos
    % see what needs to be finalized in terms of dimensionality reduction
    % idem for the remaining variables: H, ...
else
    %load(['CYG_data_raw_ind=' num2str(ind) '.mat']);
    load(['CYG_data_raw_dr_ind=' num2str(ind) '.mat']); % load only the variables needed at one point (otherwise, quite large!)
end

% figure(1),
% for i = length(uw)
%     scatter(uw{i},vw{i},'.');hold on;
% end

%% adapt to DR (see codes from Ming / ask Ming directly if necessary)
%%
ch = [1 : nChannels];
% output from the definition of the Fourier reduction scheme
H = cell(nChannels, 1);  % holographic matrices
yT = cell(nChannels, 1); % reduced data
T = cell(nChannels, 1);  % preconditioning matrix (inverse of the singular values)
Wm = cell(nChannels, 1); % masks
aW = cell(nChannels, 1); % 1./T (shall I keep this variable?)

u = cell(nChannels, 1);
v = cell(nChannels, 1);
uvidx = cell(nChannels, 1);

for i = ch % loop over the 32 channels
    
    i
    
    file_y = matfile('CYG_y.mat');
    file_pos = matfile('CYG_pos.mat');
    file_uw = matfile('CYG_uw.mat'); 
    file_vw = matfile('CYG_vw.mat');
    file_time_ = matfile('CYG_time.mat');
    file_nWw = matfile('CYG_nWw.mat');
    
    %% compute weights
    y = file_y.y(i, 1); % to be verified
    uw = file_uw.uw(i, 1);
    vw = file_vw.vw(i, 1);
    pos = file_pos.pos(i, 1);
    time_ = file_time_.time_(i, 1);
    
    % set the blocks structure
    if param_block_structure.use_manual_partitioning == 1
        param_block.size = 90000; length(uw);
        param_block.snapshot = 0;
    end
    
    u{i} = [];
    v{i} = [];
    nW{i} = [];
    uvidx{i} = [];
    
    % TO BE FIXED (only preliminary for now)
    for k = 1:numel(nSpw) % concatenate all the blocks coming from the 16 consecutive channels (loop over the 16 channels)
        % set the blocks structure (per frequency) [is aW{i} necessary?]
        param_block.pos = pos{1}{k};
        out_block = util_time_based_block_sp_ar(uw{1}{k},time_{1}{k},param_block);
        param_block_structure.partition = out_block.partition;
        aWw = util_gen_preconditioning_matrix(uw{1}{k}, vw{1}{k}, param_precond); % see if still necessary: needs to be kept in memory for DR? (no a priori)
        [u1, v1, ~, uvidx1, aW, nW1] = util_gen_block_structure(uw{1}{k}, vw{1}{k}, aWw, nWw{1}{k}, param_block_structure);
        
        % explode u1, put u1{1} concatenated in a variable, u1{2} in
        % another, ...
        
        
        % reorder the blocks along a single dimension (remove intermediate level of the frequencies)
        u{i} = [u{i}; u1];
        v{i} = [v{i}; v1];
        nW{i} = [nW{i}; nW1];
        uvidx{i} = [uvidx{i}; uvidx1];
    end
        
        % measurement operator initialization
        fprintf('Initializing the NUFFT operator\n\n');
        
        %% --- Fourier reduction --- %%
        H{i} = cell(numel(u1), 1); % one structure per block
        T{i} = cell(numel(u1), 1);
        yT{i} = cell(numel(u1), 1);
        Wm{i} = cell(numel(u1), 1);
        
        for j = 1:numel(u1) % loop over the blocks for the current channel
            [A, At, G, ~] = op_p_nufft([v1(j) u1(j)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW1(j)); 
            % this function can take cell inputs (but no need to keep multiple G matrices in memory)
            H{i}{j} = (G{1}')*G{1};
        
            % estimate threshold
            dirty2 = norm(operatorPhit(y{1}{k}(uvidx1{j}), G{1}', At) / sqrt(N));

            % no longer invoke function fourierReduction to reduce lambda functions
            % fast way of matrix probing (using psf)
            dirac2D = zeros(Ny, Nx);
            dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;
            PSF = operatorIpsf(dirac2D, A, At, H{i}{j}, [oy*Ny, ox*Nx]);
            covariancemat = FT2(PSF);
            d_mat = abs(real(covariancemat(:)));
            clear covariancemat
            %
            if param_fouRed.enable_klargestpercent
                Mask = (d_mat >= prctile(d_mat,100-param_fouRed.klargestpercent));
            elseif param_fouRed.enable_estimatethreshold
                % Embed the noise
                noise1 = param_fouRed.sigma_noise * (randn(size(G{1}, 1),1) + 1j * randn(size(G{1}, 1), 1));
                rn = FT2(At(G{1}'*noise1));  % Apply F Phi
                th_dirty = param_fouRed.gamma * std(rn(:)) / dirty2;
                fprintf('\nThe estimate threshold using ground truth is %e \n', th);
                Mask = (d_mat >= th_dirty);
            end
            d_mat = d_mat(Mask);
            %
            T{i}{j} = max(param_fouRed.diagthresholdepsilon, d_mat);  % This ensures that inverting the values will not explode in computation
            T{i}{j} = 1./sqrt(T{i}{j});
            Wm{i}{j} = Mask;

            % reduce the data block
            yT{i}{j} = dataReduce(y{1}{k}(uvidx1{j}), G{1}', At, T{i}{j}, Wm{i}{j});
        end
    end
    % missing: definition of epsilon... [to be checked with Ming]
    %% ---
    
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
