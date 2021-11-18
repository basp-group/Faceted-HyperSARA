clc; clear all; close all
format compact;

addpath lib/faceted-wavelet-transform/src

%%
nchannels = 1;
nblocks_per_channel = 2;
block_size = 5;

% this needs to be revised / generalized
nfacets = 0;
ncores_channel = nchannels;
ncores_block = nblocks_per_channel;
ncores_data = ncores_channel * ncores_block;
numworkers = ncores_data + nfacets;

%%
cirrus_cluster = parcluster('local');
cirrus_cluster.NumWorkers = numworkers;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
parpool(cirrus_cluster, numworkers);


%%
% rg_c = split_range(ncores_channel, nchannels);

spmd
    % on any data worker
    if labindex > nfacets
        frequency_id = 1 + floor((labindex - nfacets - 1)/ncores_block);
        
        if mod(labindex - nfacets - 1, ncores_block) == 0
            % on root frequency data worker
            y = labindex * ones(nblocks_per_channel*block_size, 1);

            % send portions of y to each data block worker (consecutive labindex)
            % first portion will be handled by the frequency root itself
            local_y = y(1:block_size);
            for k = 2:ncores_block
                labSend(y((k-1)*block_size + 1 : k*block_size), labindex + (k-1))
            end       
        else
            % on data block workers
            root_frequency_worker = nfacets + 1 + (frequency_id - 1) * ncores_block + mod(labindex - nfacets, ncores_block);
            local_y = labReceive(root_frequency_worker);
        end
    end
end

%% root data worker (in charge of the fft)

