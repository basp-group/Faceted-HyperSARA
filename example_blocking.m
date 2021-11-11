clc; clear all; close all
format compact;

addpath lib/faceted-wavelet-transform/src

%% base characteristics
nfacets = 1;
nchannels = 3;

% channel groups
nchannel_groups = 2; % > nchannels
ncores_per_channel_group = [2, 1]; % (array of size [nchannel_groups, 1])

% data blocks
nblocks_per_channel = [1, 2, 1]; % (array of size [nchannels, 1])

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
        
        % list of frequencies handled locally (return empty if none if handled)
        % split evenly between workers in the group
        local_channels = local_split_range(ncores_per_channel_group(frequency_group_id), frequency_group_nchannels, frequency_subgroup_id);
        global_channels = global_start_frequency + local_channels - 1;
        
        % total number of data blocks handled in the frequency group
        frequency_group_nblocks = sum(nblocks_per_channel(channels_groups(frequency_group_id,1):channels_groups(frequency_group_id,2))); 
        % = nblocks_per_frequency_group(frequency_group_id)
        
        % identify global id of the channels associated to each block handled
        % locally
        
        % get id of the workers responsible for those channnels
        
        % identify global block id 
    end
end

