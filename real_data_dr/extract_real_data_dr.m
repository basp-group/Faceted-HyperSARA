%% Real data extraction
addpath ../lib/utils/
addpath ../fouRed/
addpath ../../real_data_dr
addpath ../lib/operators
addpath ../lib/nufft

visibility_file_name = 'CYG_data_raw_ind=';
param_real_data.image_size_Nx = 2560;
param_real_data.image_size_Ny = 1536;
nSpw = 2;           % number of spectral channels per MS file
nChannels = 2*nSpw; % total number of "virtual" channels (i.e., after 
                    % concatenation) for the real dataset considered
nBlocks = 9;        % number of data blocks (needs to be known beforehand,
                    % quite restrictive here), change l.70 accordingly
klargestpercent = 20;
extract_raw_data = true;
% param_real_data.pixel_size = 0.3; % in arcsec
FT2 = @(x) fftshift(fft2(ifftshift(x)));

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

%% Fourier reduction parameters
param_fouRed.enable_klargestpercent = 1;
param_fouRed.klargestpercent = klargestpercent;
param_fouRed.enable_estimatethreshold = 0;
param_fouRed.gamma = 3;             % By using threshold estimation, the optimal theshold reads as gamma * sigma / ||x||_2
param_fouRed.diagthresholdepsilon = 1e-10;
param_fouRed.covmatfileexists = 0;
param_fouRed.covmatfile = 'covariancemat.mat';
param_fouRed.fastCov = 1;

%% Block structure
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

if param_block_structure.use_manual_partitioning == 1
    param_block.size = 90000; % to be changed (make sure nBlocks in the end...)
    param_block.snapshot = 0;
end

%% Load/extract real data from band-interleaved .mat files
if extract_raw_data
    % visibilities
    y = cell(nChannels, 1);
    for l = 1:numel(y)
        y{l} = cell(nBlocks, 1);
    end
    save('CYG_y.mat', 'y', '-v7.3');
    clear y;
        
    % u
    u = cell(nChannels, 1);
    for l = 1:numel(u)
        u{l} = cell(nBlocks, 1);
    end
    save('CYG_u.mat', 'u', '-v7.3');
    clear u;
    
    % v
    v = cell(nChannels, 1);
    for l = 1:numel(v)
        v{l} = cell(nBlocks, 1);
    end
    save('CYG_v.mat', 'v', '-v7.3');
    clear v;
    
    % nWw: caling NUFFT
    nW = cell(nChannels, 1);
    for l = 1:numel(nW)
        nW{l} = cell(nBlocks, 1);
    end
    save('CYG_nW.mat', 'nW', '-v7.3');
    clear nW;
    
%     % position of data for different configs inside each visibility vector
%     % (e.g., A, B, C ...)
%     pos = cell(nChannels, 1);
%     for l = 1:numel(pos)
%         pos{l} = cell(nBlocks, 1);
%     end
%     save('CYG_pos.mat', 'pos', '-v7.3');
%     clear pos;    
    
%     % acquisition time
%     time = cell(nChannels, 1);
%     for l = 1:numel(time)
%         time{l} = cell(nBlocks, 1);
%     end
%     save('CYG_time.mat', 'time', '-v7.3');
%     clear time;
    
    new_file_y = matfile('CYG_y.mat', 'Writable', true);
    new_file_u = matfile('CYG_u.mat', 'Writable', true);
    new_file_v = matfile('CYG_v.mat', 'Writable', true);
    new_file_nW = matfile('CYG_nW.mat', 'Writable', true);
    
    for l = 1:2*nSpw
        for n = 1:nSpw
            filename = [visibility_file_name, num2str(n), '.mat'];
            file = matfile(filename);
            y = cell2mat(file.y(1,l));
            pos = cell2mat(file.pos(1,l));
            time = cell2mat(file.time_(1,l));
            uw = cell2mat(file.uw(1,l));
            vw = cell2mat(file.vw(1,l));
            nWw = cell2mat(file.nWw(1,l));
            
            % blocking
            % set the blocks structure (per frequency)
            param_block.pos = pos;
            out_block = util_time_based_block_sp_ar(uw, time,param_block);
            param_block_structure.partition = out_block.partition;
            aWw = util_gen_preconditioning_matrix(uw, vw, param_precond);
            [u1, v1, ~, uvidx1, aW, nW1] = util_gen_block_structure(uw, vw, aWw, nWw, param_block_structure);
            
            % concatenate the visibility frequencies and the associated u/v
            % points (make sure same number of blocks as )
            u_tmp = new_file_u.u(l,1);
            v_tmp = new_file_v.v(l,1);
            y_tmp = new_file_y.y(l,1);
            nW_tmp = new_file_nW.nW(l,1);
            for m = 1:numel(u1)
                u_tmp{1}{m} = [u_tmp{1}{m}; u1{m}];
                v_tmp{1}{m} = [v_tmp{1}{m}; v1{m}];
                y_tmp{1}{m} = [y_tmp{1}{m}; y(uvidx1{m})];
                nW_tmp{1}{m} = [nW_tmp{1}{m}; nW1{m}];
            end
            new_file_u.u(l,1) = u_tmp;
            new_file_v.v(l,1) = v_tmp;
            new_file_y.y(l,1) = y_tmp;
            new_file_nW.nW(l,1) = nW_tmp;
        end
    end
else
    new_file_y = matfile('CYG_y.mat');
    new_file_u = matfile('CYG_u.mat');
    new_file_v = matfile('CYG_v.mat');
    new_file_nW = matfile('CYG_nW.mat');
end

%% Estimate epsilon with NNLS on each data block
... to be implemented here ...

%% Define DR operators / reduce data blocks
%%
% output from the definition of the Fourier reduction scheme
H = cell(nChannels, 1);  % holographic matrices
yT = cell(nChannels, 1); % reduced data
T = cell(nChannels, 1);  % preconditioning matrix (inverse of singular values)
Wm = cell(nChannels, 1); % masks for the Fourier plane corresponding to 
                         % data blocks

% loop over virtual channels
for l = 1:nChannels
    
    % measurement operator initialization
    fprintf('Initializing the NUFFT operator\n\n');
    
    %% --- Fourier reduction --- %%  
    y_tmp = new_file_y.y(l,1);
    u_tmp = new_file_u.u(l,1);
    v_tmp = new_file_v.v(l,1);
    nW_tmp = new_file_nW.nW(l,1);
    
    H{l} = cell(numel(u_tmp{1}), 1);
    T{l} = cell(numel(u_tmp{1}), 1);
    yT{l} = cell(numel(u_tmp{1}), 1);
    Wm{l} = cell(numel(u_tmp{1}), 1);
    
    % loop over the data blocks within each "virtual" channel
    for j = 1:numel(u_tmp{1})
        
        [A, At, G, ~] = op_p_nufft([v_tmp{1}(j) u_tmp{1}(j)], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2], nW_tmp{1}(j));
        % note: the above function can take cell inputs (but no need to keep multiple G matrices in memory)
        H{l}{j} = (G{1}')*G{1}; % pprogressively write to disk? (huge...)
        
        % estimate threshold
        dirty2 = norm(operatorPhit(y_tmp{1}{j}, G{1}', At) / sqrt(N));
        
        % fast matrix probing (using psf)
        dirac2D = zeros(Ny, Nx);
        dirac2D(ceil((Ny+1)/2), ceil((Nx+1)/2)) = 1;
        PSF = operatorIpsf(dirac2D, A, At, H{l}{j}, [oy*Ny, ox*Nx]);
        covariancemat = FT2(PSF);
        d_mat = abs(real(covariancemat(:)));
        clear covariancemat
        %
        if param_fouRed.enable_klargestpercent
            Mask = (d_mat >= prctile(d_mat,100-param_fouRed.klargestpercent));
        elseif param_fouRed.enable_estimatethreshold
            % embed the noise
            noise1 = param_fouRed.sigma_noise * (randn(size(G{1}, 1),1) + 1j * randn(size(G{1}, 1), 1));
            rn = FT2(At(G{1}'*noise1));  % apply F Phi
            th_dirty = param_fouRed.gamma * std(rn(:)) / dirty2;
            fprintf('\nThe estimate threshold using ground truth is %e \n', th);
            Mask = (d_mat >= th_dirty);
        end
        d_mat = d_mat(Mask);
        %
        T{l}{j} = max(param_fouRed.diagthresholdepsilon, d_mat);  % ensures that inverting the values will not explode in computation
        T{l}{j} = 1./sqrt(T{l}{j});
        Wm{l}{j} = Mask;
        
        % reduce the data block
        yT{l}{j} = dataReduce(v_tmp{1}{j}, G{1}', At, T{l}{j}, Wm{l}{j});
    end
    % missing: definition of epsilon... [to be checked with Ming]: do NNLS
    % beforehand...
    %% ---
end

% %% Save data
% if save_data
%     save('CYG_data.mat','-v7.3', 'G', 'W', 'aW', 'yb');
%     save('CYG_y.mat','-v7.3', 'y');
% end
%
% %%
% if save_full_operator && exist('Gw','var')
%     save('CYG_Gw.mat','-v7.3', 'Gw');
% end
%
% %% Free memory
% if free_memory
%     clear y u v uv_mat uv_mat1 uvidx uw vw ant1 ant2 aWw nW nWw out_block;
% end
