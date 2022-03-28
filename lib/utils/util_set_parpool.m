function hpc_cluster = util_set_parpool(solver, ncores_data, Q, parcluster_name)
    % function to launch parpool before the solver.
    %
    % Parameters
    % ----------
    % algo_version : [string]
    %      Solver used, 'fhs', 'hs', or 'sara'
    % ncores_data : [int]
    %     number of data workers
    % Q : [int]
    %     number of facet workers
    % parcluster_name : [string]
    %     name of the parcluster profile to use, default "local"
    %
    % Returns
    % -------
    % hpc_cluster
    %     [description]

    switch solver
        case 'sara'
            numworkers = ncores_data;
        case 'hs'
            % 2 facets workers (main session), ncores_data data workers
            numworkers = ncores_data;
        case 'fhs'
            % Q facets workers, ncores_data: data workers
            numworkers = Q + ncores_data;
    end

    hpc_cluster = parcluster(parcluster_name); 
    hpc_cluster.NumWorkers = numworkers;
    hpc_cluster.NumThreads = 1;
  
    % create temp. dir. if slurm parcluster profile is used
    if ~strcmp(parcluster_name,'local')
        
        tempdir = [pwd,filesep,'tmp_files',filesep];
        mkdir(tempdir)
        try tempdir=[tempdir, getenv('SLURM_JOB_ID')];
        end
        mkdir(tempdir)
        hpc_cluster.JobStorageLocation = tempdir;
    end
   
    % % start the matlabpool with the requested workers
    parpool(hpc_cluster, numworkers);
    
    % wavelet 
    dwtmode('zpd', 'nodisp');
    spmd
        dwtmode('zpd', 'nodisp');
    end

end
