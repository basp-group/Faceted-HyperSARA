clc; clear all; close all
format compact;

addpath lib/faceted-wavelet-transform/src

%% Base characteristics (to be provided by the user)
nfacets = 0;
nchannels = 3;
ncores_data = 2; % 4

% data blocks
nblocks_per_channel = [2, 3, 1]; % (array of size [nchannels, 1])
if numel(nblocks_per_channel) < nchannels
    error('nblocks_per_channel needs to contain nchannels elements')
end

% total number of cores required
numworkers = ncores_data + nfacets;

%%

% list of channels handled by each worker
if ncores_data < nchannels
    channels_groups = split_range(ncores_data, nchannels);
else
    channels_groups = split_range(nchannels, nchannels);
    % ! 0 as the stop frequency means that the worker is not responsible for
    % any given frequency
    channels_groups = [channels_groups; [nchannels*ones(ncores_data - nchannels, 1), zeros(ncores_data - nchannels, 1)]];
end

% list of blocks (unrolled indices) handled by each worker
% ! the case ncores_data > nblocks should never occur (does not make amy sense)
nblocks = sum(nblocks_per_channel );
if ncores_data > nblocks
    error('The number of data cores, ncores_data, should be greater than or equal to the total number of blocks');
end
block_groups = split_range(ncores_data, nblocks);

% channel id of each block
% TODO: see if computation can be simplified
% TODO: see if the info can be computed locally, w/o any interaction with
% other workers
c1 = cumsum(nblocks_per_channel);
c0 = [0, c1];
blocks_channel_id = zeros(nblocks, 1);
blocks_channel_worker_id = zeros(nblocks, 1);
blocks_within_channel_id = zeros(nblocks, 1);
for b = 1:nblocks
    blocks_channel_id(b) = find(c1 >= b, 1, 'first');
    blocks_channel_worker_id(b) = find(channels_groups(:, 2) >= blocks_channel_id(b), 1, 'first'); % id of the worker responsible for the channel
    blocks_within_channel_id(b) = b - c0(blocks_channel_id(b));
end

%%
y = cell(nchannels, 1);
gt_y = zeros(nchannels, 1);

for l = 1:nchannels
    y{l} = cell(nblocks_per_channel(l), 1);
    for b = 1:nblocks_per_channel(l)
        y{l}{b} = rand(5, 1);
        gt_y(l) = gt_y(l) + norm(y{l}{b})^2;
    end
    
end

%% parallelization
cirrus_cluster = parcluster('local');
cirrus_cluster.NumWorkers = numworkers;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
parpool(cirrus_cluster, numworkers);

%%

% TODO: group communications of the Fourier plane (i.e., if multiple channels involved)
spmd
    % on any data worker
    if labindex > nfacets
        
        % id among all data workers
        data_id = labindex - nfacets;
        
        % list of channels handled on current worker
        if data_id > ncores_data
            local_channels = [nchannels, 0];
            local_nchannels = 0;
        else
            local_channels = local_split_range(ncores_data, nchannels, data_id);
            local_nchannels = local_channels(2) - local_channels(1) + 1;
        end

        % list of blocks handled on current worker
        local_blocks = local_split_range(ncores_data, nblocks, data_id);
        local_nblocks = local_blocks(2) - local_blocks(1) + 1;
        
        % identify channel associated to each local block (retrieved from central node...)
        local_blocks_channel_id = blocks_channel_id(local_blocks(1):local_blocks(2));
        
        local_blocks_within_channel_id = blocks_within_channel_id(local_blocks(1):local_blocks(2));
        
        % for each local block, identify on which worker the frequency it
        % derives from is handled
        local_blocks_channel_worker_id = blocks_channel_worker_id(local_blocks(1):local_blocks(2));
        
        % identify number of distinct channels to which local blocks are
        % associated
        local_distinct_channels_of_blocks = numel(unique(local_blocks_channel_id));
        local_y = zeros(local_distinct_channels_of_blocks, 1);
        
        % check consistency 
        l = 1;
        m = min(local_blocks_channel_id);
        for b = 1:local_nblocks
            l = local_blocks_channel_id(b) - m + 1;
            local_y(l) = local_y(l) + norm(y{local_blocks_channel_id(b)}{local_blocks_within_channel_id(b)})^2;
        end
        
        % identify values which need to be communicated across workers
        % (i.e., whenever blocks_channel_worker_id(b) differs from data_id)
        
%         for l = 1:local_nchannels
%             %local_y{l} = cell(local_nblocks, 1);
%             for b = 1:local_nblocks
%                 %local_y{l}{b} = y{local_channels + l - 1}{local_blocks_within_channel_id{b}};
%                 % blocks from outer channels (i.e., stored on other workers)
%                 if local_blocks_channel_id(b) == l
%                     local_y(l) = local_y(l) + norm(y{local_blocks_channel_id(b)}{local_blocks_within_channel_id(b)})^2;
%                 end
%             end
%         end

        % send appropriate value of the sum
        
        % receive values (if needed)
        % ! receive mode only if blocks from the channels handled locally
        % have been stored elsewhere
    
    end
end

