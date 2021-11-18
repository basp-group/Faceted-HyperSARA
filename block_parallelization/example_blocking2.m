clc; clear all; close all
format compact;

addpath lib/faceted-wavelet-transform/src

%% base characteristics
nfacets = 0;
nchannels = 3;

% channel groups
nchannel_groups = 1; % > nchannels
ncores_per_channel_group = [4]; % (array of size [nchannel_groups, 1])

% data blocks
nblocks_per_channel = [2, 2, 1]; % (array of size [nchannels, 1])

% total number of cores required
% this needs to be revised / generalized
ncores_data = sum(ncores_per_channel_group);
numworkers = ncores_data + nfacets;

%% create mock data
y = cell(nchannels, 1);

for l = 1:nchannels
    y{l} = cell(nblocks_per_channel(l), 1);
    for b = 1:nblocks_per_channel(l)
        y{l}{b} = rand(5, 1);
    end
end

%%
cirrus_cluster = parcluster('local');
cirrus_cluster.NumWorkers = numworkers;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
parpool(cirrus_cluster, numworkers);

%%

channels_groups = split_range(nchannel_groups, nchannels);
nblocks_per_frequency_group = zeros(nchannel_groups, 1);
for l = 1:nchannel_groups
    nblocks_per_frequency_group(l) = sum(nblocks_per_channel(channels_groups(l, 1) : channels_groups(l, 2)));
end
c0 = cumsum(nblocks_per_frequency_group);
c1 = cumsum(nblocks_per_channel);

%%

spmd
    % on any data worker
    if labindex > nfacets
        
        % id among all data workers
        data_id = labindex - nfacets;
        % frequency group
        frequency_group_id = find(cumsum(ncores_per_channel_group) >= data_id, 1, 'first'); % -1
        global_start_frequency = channels_groups(frequency_group_id, 1);
        % total number of channels in the frequency group
        frequency_group_nchannels = channels_groups(frequency_group_id, 2) - channels_groups(frequency_group_id, 1) + 1; 
        % frequency subgroup index (id among workers in the same frequency group)
        frequency_subgroup_id = data_id - sum(ncores_per_channel_group(1:frequency_group_id - 1));
        % -> ideally, should be 0 for cores handling only blocks from a channels ? (and 1 for the root of the frequency group)
        % frequency root workers are such that frequency_subgroup_id == 1
        
        % list of frequencies handled locally (return empty if none is handled)
        % split evenly between workers in the group
        if frequency_group_nchannels >= ncores_per_channel_group(frequency_group_id)
            % more channel than cores for the group
            local_channels = local_split_range(ncores_per_channel_group(frequency_group_id), frequency_group_nchannels, frequency_subgroup_id);
            global_channels = global_start_frequency + local_channels - 1;
            local_nchannels = local_channels(2) - local_channels(1) + 1;
        else
            % TODO: to be checked against all possible configurations
            % fewer channels than cores for the group
            % create channel identifiers only for workers fully responsible
            % for a channel
            if data_id <= frequency_group_nchannels
                local_channels = local_split_range(frequency_group_nchannels, frequency_group_nchannels, frequency_subgroup_id);
                global_channels = global_start_frequency + local_channels - 1;
                local_nchannels = local_channels(2) - local_channels(1) + 1;
            else
                local_channels = 0; % [1, frequency_group_nchannels]; % ! be careful here
                global_channels = 0; % global_start_frequency; % ! idem, to be tested
                local_nchannels = 0;
            end
        end
        
        % total number of data blocks handled in the frequency group
        nblocks_per_channel_frequency_group = nblocks_per_channel(channels_groups(frequency_group_id,1):channels_groups(frequency_group_id,2));
        frequency_group_nblocks = sum(nblocks_per_channel_frequency_group); % = nblocks_per_frequency_group(frequency_group_id)
        
        % unrolled index of the blocks to be handled on the current worker
        % ! there will never be more cores than blocks in a channel group
        % (this possibility needs to be ruled out, and trigger an error early)
        local_unrolled_block_ids = local_split_range(ncores_per_channel_group(frequency_group_id), frequency_group_nblocks, frequency_subgroup_id);
        local_nblocks = local_unrolled_block_ids(2)-local_unrolled_block_ids(1)+1; % number of blocks handled on current worker
        global_unrolled_block_ids = local_unrolled_block_ids + sum(nblocks_per_frequency_group(1:frequency_group_id-1));
        
        % associated frequency for each block (id within frequency subgroup)
        % ! problem here
%         c = cumsum(nblocks_per_channel_frequency_group);
%         blocks_frequency_group_id = zeros(local_nblocks, 1);
%         blocks_id_within_frequency = zeros(local_nblocks, 1); % [local_nblocks, 1]
%         for b = 1:local_nblocks
%             blocks_frequency_group_id(b) = find(c >= data_id - c0(frequency_group_id) - local_unrolled_block_ids(1) + b - 1, 1, 'first');
%             blocks_id_within_frequency(b) = data_id - sum(c(1:blocks_frequency_group_id(b)-1));
%         end

        c = cumsum(nblocks_per_channel_frequency_group);
        blocks_global_frequency_id = zeros(local_nblocks, 1); % -> retrieve id of worker to communicate with
        blocks_local_frequency_id = zeros(local_nblocks, 1);
        blocks_within_frequency_id = zeros(local_nblocks, 1); % [local_nblocks, 1]
        for b = 1:local_nblocks
            blocks_global_frequency_id(b) = find(c1 >= global_unrolled_block_ids(1) + b - 1, 1, 'first'); % ok
            
            % problems here
            %blocks_channel_root_worker(b) = 
            %blocks_channel_root_local_frequency_id(b) = % local frequency id on the core responsible for the channels
            %blocks_within_frequency_id(b) = find(c >= blocks_channel_root_local_frequency_id(b), 1, 'first'); % problem here for now
            
            %blocks_local_frequency_id(b) = blocks_global_frequency_id(b) - (global_start_frequency + b - 1) + 1; % problem here
            %blocks_within_frequency_id(b) = find(c >= blocks_local_frequency_id(b), 1, 'first'); % problem here
        end
        
        % get id of the workers responsible for those channnels associated
        % with each block handled locally (in case blocks stored on a worker different from the
        % one responsible for the associated channel)
    end
end

