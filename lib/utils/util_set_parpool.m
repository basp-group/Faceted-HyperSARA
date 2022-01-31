function cirrus_cluster = util_set_parpool(algo_version, ncores_data, Q, flag_cirrus)
    % [extended_summary]
    %
    % Parameters
    % ----------
    % algo_version : [type]
    %     [description]
    % ncores_data : [type]
    %     [description]
    % Q : [type]
    %     [description]
    % flag_cirrus : [type]
    %     [description]
    %
    % Returns
    % -------
    % [type]
    %     [description]

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

    cirrus_cluster = parcluster('local');
    cirrus_cluster.NumWorkers = numworkers;
    cirrus_cluster.NumThreads = 1;
    ncores = cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads;
    if cirrus_cluster.NumWorkers * cirrus_cluster.NumThreads > ncores
        exit(1);
    end
    % explicitly set the JobStorageLocation to the temp directory that was created in your sbatch script
    if flag_cirrus
        cirrus_cluster.JobStorageLocation = strcat('/lustre/home/sc004/', getenv('USER'), '/', getenv('SLURM_JOB_ID'));
    end
    % maxNumCompThreads(param.num_workers);
    parpool(cirrus_cluster, numworkers);
    % % start the matlabpool with maximum available workers
    % % control how many workers by setting ntasks in your sbatch script
    % parpool(cirrus_cluster, str2num(getenv('SLURM_CPUS_ON_NODE')))
    dwtmode('zpd', 'nodisp');
    spmd
        dwtmode('zpd', 'nodisp');
    end

end
