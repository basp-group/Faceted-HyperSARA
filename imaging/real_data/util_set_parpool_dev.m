function hpc_cluster = util_set_parpool_dev(algo_version, ncores_data, Q, flag_hpc)

switch algo_version
    case 'sara'
        numworkers = ncores_data;
    case 'hs'
        % 2 facets workers (main session), ncores_data data workers
        numworkers = ncores_data;
    case 'fhs'
        % Q facets workers, ncores_data: data workers
        numworkers = Q + ncores_data;
end

try delete(gcp('nocreate'));
end

if flag_hpc
       fprintf('\nINFO: will open a parpool of size %d  using a slurm parcluster',numworkers)
       hpc_cluster = parcluster('mySlurmProfileSingleThread'); % update the parcluster name
       hpc_cluster.ResourceTemplate = '--export=ALL --ntasks=^N^ --cpus-per-task=^T^ --tasks-per-node=36 --account=ec110-guest  --time=06:00:00 --partition=standard --qos=standard';
else,  hpc_cluster = parcluster('local');
end

hpc_cluster.NumWorkers = numworkers;
hpc_cluster.NumThreads = 1;
ncores = hpc_cluster.NumWorkers * hpc_cluster.NumThreads;
if hpc_cluster.NumWorkers * hpc_cluster.NumThreads > ncores
    exit(1);
end
% explicitly set the JobStorageLocation to the temp directory that was created in your sbatch script
if flag_hpc
    pathtemp = [pwd,filesep,'tmp_parpool',filesep];
    mkdir(pathtemp);
    pathtemp = [pathtemp, getenv('SLURM_JOB_ID')];
    mkdir(pathtemp);
    hpc_cluster.JobStorageLocation = pathtemp;
end


parpool(hpc_cluster, numworkers);

% % start the matlabpool with maximum available workers
% % control how many workers by setting ntasks in your sbatch script
% parpool(hpc_cluster, str2num(getenv('SLURM_CPUS_ON_NODE')))
dwtmode('zpd', 'nodisp');
spmd
    dwtmode('zpd', 'nodisp');
end

end
