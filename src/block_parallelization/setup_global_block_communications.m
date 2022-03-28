function [channels_groups, block_groups, blocks_channel_worker_id, ...
    blocks_within_channel_id] = ...
    setup_global_block_communications(nchannels, nblocks_per_channel, ...
    ncores_data, channel_core, block_core)
%% basic checkup
if numel(nblocks_per_channel) < nchannels
    error('nblocks_per_channel needs to contain nchannels elements');
end

if numel(channel_core) < nchannels
    error('channel_core needs to contain nchannels elements');
end

% list of blocks (unrolled indices) handled by each worker
% ! the case ncores_data > nblocks should never occur (does not make any sense)
nblocks = sum(nblocks_per_channel);
if ncores_data > nblocks
    error('The number of data cores, ncores_data, should be greater than or equal to the total number of blocks');
end

if numel(block_core) < nblocks
    error('block_core needs to contain sum(nblocks_per_channel) elements');
end

%%
% list of channels handled by each worker
channels_groups = zeros(ncores_data, 2);
% TODO: improve this part
for c = 1:min(ncores_data, nchannels)
    a = find(channel_core == c, 1, 'first');
    if ~isempty(a)
        channels_groups(c, 1) = a;
    else
        channels_groups(c, 1) = 1;
    end
    a = find(channel_core == c, 1, 'last');
    if ~isempty(a)
        channels_groups(c, 2) = a;
    end
end

if nchannels < ncores_data
    % ! 0 as the stop frequency means that the worker is not responsible for
    % any given frequency
    channels_groups(nchannels + 1:end, 1) = nchannels;
end

% TODO: improve this part
block_groups = zeros(ncores_data, 2);
for c = 1:min(ncores_data, nchannels)
    a = find(block_core == c, 1, 'first');
    if ~isempty(a)
        block_groups(c, 1) = a;
    else
        block_groups(c, 1) = 1;
    end
    a = find(block_core == c, 1, 'last');
    if ~isempty(a)
        block_groups(c, 2) = a;
    end
end

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

end
