function cirrus_cluster = util_set_parpool_dev_ad(algo_version, ncores_data, Q, flag_cirrus)
  
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


delete(gcp('nocreate'));

if flag_cirrus
       cirrus_cluster =parcluster('mySlurmProfileSingleThread');  %parcluster('local');
       cirrus_cluster.ResourceTemplate = '--export=ALL --ntasks=^N^ --cpus-per-task=^T^ --tasks-per-node=36 --account=ec110-guest --time=96:00:00 --partition=standard --qos=standard';
else
       cirrus_cluster =parcluster('local');
end

cirrus_cluster.NumWorkers = numworkers;
cirrus_cluster.NumThreads = 1;
ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
if cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads > ncores
    exit(1);
end
% explicitly set the JobStorageLocation to the temp directory that was created in your sbatch script
%if flag_cirrus
    path=[strcat('/lustre/home/sc004/', getenv('USER'),'/matlab_jobs_tmp/')];
    mkdir(path)
    path_dummy = [path, getenv('SLURM_JOB_ID')];
    mkdir(path_dummy)
    cirrus_cluster.JobStorageLocation = path_dummy;% strcat('/lustre/home/sc004/', getenv('USER'),'/', getenv('SLURM_JOB_ID'));
%end


% maxNumCompThreads(param.num_workers);

parpool(cirrus_cluster, numworkers);


% % start the matlabpool with maximum available workers
% % control how many workers by setting ntasks in your sbatch script
% parpool(cirrus_cluster, str2num(getenv('SLURM_CPUS_ON_NODE')))
dwtmode('zpd','nodisp')
spmd
    dwtmode('zpd','nodisp')
end

end
