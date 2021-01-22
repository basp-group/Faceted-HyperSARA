function [aruntime, vruntime, acpu_time, vcpu_time, total_runtime, total_cpu_time, iteration_number] = get_timing_sara(results_filename, metric_filename)
    % results_file : dictionary a strings, pointing to the mat corresponding to
    % the different (spectral faceting) sub-cubes.
    
    Qc = numel(results_filename); % number of spectral cubes (channels for SARA)
    
    sum_runtime_sqr = 0;
    sum_l11_sqr = 0;
    sum_master_sqr = 0;
    sum_data_sqr = 0;
    sum_cpu_sqr = 0;
    
    total_runtime = 0;  % assuming all sub-cubes are running in parallel
    total_cpu_time = 0; % 
    aruntime = 0;       % average time (per iteration)
    atime_l11 = 0;    
    atime_master = 0;
    atime_data = 0;     
    iteration_number = 0;
    ncores_l11 = 9; % number of SARA dictionaries
    
    for ind = 1:Qc
        load(results_filename{ind}, 't_l11', 't_master', 't_data', 'end_iter')
        
        % compute mean and variance over all the files? just 1 for now
        aruntime = aruntime + sum(end_iter(end_iter > 0)); % average runtime per iteration, over all sub-problems
        sum_runtime_sqr = sum_runtime_sqr + sum(end_iter(end_iter > 0).^2);
    
        atime_l11 = atime_l11 + sum(t_l11(t_l11 > 0));
        atime_master = atime_master + sum(t_master(t_master > 0));
        atime_data = atime_data + sum(t_data(t_data > 0));
        sum_l11_sqr = sum_l11_sqr + (ncores_l11^2)*sum(t_l11(t_l11 > 0).^2);
        sum_master_sqr = sum_master_sqr + sum(t_master(t_master > 0).^2);
        sum_data_sqr = sum_data_sqr + sum(t_data(t_data > 0).^2);
    
        % average number of iterations over all sub-problems
        iteration_number = iteration_number + sum(end_iter > 0);
        total_runtime = max(total_runtime, sum(end_iter(end_iter > 0)));
        total_cpu_time = total_cpu_time + ncores_l11*sum(t_l11(t_l11 > 0)) ...
            + sum(t_master(t_master > 0)) + sum(t_data(t_data > 0));
        sum_cpu_sqr = sum_cpu_sqr + sum((ncores_l11*t_l11(t_l11 > 0) ...
            + t_master(t_master > 0) + t_data(t_data > 0)).^2);
    end
    
    % iteration_number = round(iteration_number / Qc);
    % only report average iteration number over all sub-problems
    % average per iteration per sub-cube
    aruntime = aruntime/iteration_number;
    atime_l11 = atime_l11/iteration_number;
    atime_master = atime_master/iteration_number;
    atime_data = atime_data/iteration_number;
    %
    vruntime = (sum_runtime_sqr - iteration_number*aruntime^2)/(iteration_number - 1);
    vtime_l11 = (sum_l11_sqr - iteration_number*atime_l11^2)/(iteration_number - 1);
    vtime_master = (sum_l11_sqr - iteration_number*atime_master^2)/(iteration_number - 1);
    vtime_data = (sum_data_sqr - iteration_number*atime_data^2)/(iteration_number - 1);
    %
    acpu_time = total_cpu_time/iteration_number;
    vcpu_time = (sum_cpu_sqr - iteration_number*acpu_time^2)/(iteration_number - 1);
    
    for k = 1:numel(Qc)
        fprintf("Qc = %i, iteration_number = %i \n", ...
            Qc, iteration_number/Qc)
        fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
            aruntime, sqrt(vruntime), acpu_time, sqrt(vcpu_time));
        fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
            total_runtime/3600, total_cpu_time/3600)
    end
    
    %% Saving results
    save(metric_filename, '-v7.3', 'aruntime', 'vruntime', 'acpu_time', ...
        'vcpu_time', 'iteration_number', ...
        'atime_l11', 'atime_master', 'atime_data', ...
        'vtime_l11', 'vtime_master', 'vtime_data', ...
        'total_runtime', 'total_cpu_time');
    end
    