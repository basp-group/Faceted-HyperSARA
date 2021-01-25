function [aruntime, vruntime, acpu_time, vcpu_time, total_runtime, total_cpu_time, iteration_number, astd_res, akurt_res] = get_metric_facetedHypersara(results_filename, metric_filename, ncores_facets)
% results_file : dictionary a strings, pointing to the mat corresponding to
% the different (spectral faceting) sub-cubes.

Qc = numel(results_filename); % number of spectral cubes

sum_runtime_sqr = 0;
sum_facet_sqr = 0;
sum_data_sqr = 0;
sum_cpu_sqr = 0;

total_runtime = 0;  % assuming all sub-cubes are running in parallel
total_cpu_time = 0; % 
aruntime = 0;       % average time (per iteration)
atime_facet = 0;    
atime_data = 0;     
iteration_number = 0;

std_res = [];
kurt_res = [];
skew_res = [];

for ind = 1:Qc
    load(results_filename{ind}, 't_facet', 't_data', 'end_iter', 'param', 'res')
    ncores_data = param.num_workers - ncores_facets; % number of cpu assigned to the data fidelity terms
    
    std_res = [std_res; std(res,0,[1,2])];
    kurt_res = [kurt_res; kurtosis(res,0,[1,2])];
    skew_res = [skew_res;skewness(res,0,[1,2])];
    
    % compute mean and variance over all the files? just 1 for now
    aruntime = aruntime + sum(end_iter(end_iter > 0)); % average runtime per iteration, over all sub-problems
    sum_runtime_sqr = sum_runtime_sqr + sum(end_iter(end_iter > 0).^2);

    atime_facet = atime_facet + sum(t_facet(t_facet > 0));
    atime_data = atime_data + sum(t_data(t_data > 0));
    sum_facet_sqr = sum_facet_sqr + (ncores_facets^2)*sum(t_facet(t_facet > 0).^2);
    sum_data_sqr = sum_data_sqr + (ncores_data^2)*sum(t_data(t_data > 0).^2);

    % average number of iterations over all sub-problems
    iteration_number = iteration_number + sum(end_iter > 0);
    total_runtime = max(total_runtime, sum(end_iter(end_iter > 0)));
    total_cpu_time = total_cpu_time + ncores_facets*sum(t_facet(t_facet > 0)) ...
        + ncores_data*sum(t_data(t_data > 0));
    sum_cpu_sqr = sum_cpu_sqr + sum((ncores_facets*t_facet(t_facet > 0) ...
        + ncores_data*t_data(t_data > 0)).^2);
end

% average std / kurtosis for the residual
astd_res = mean(std_res);
akurt_res = mean(kurt_res);
askew_res = mean(skew_res);

% iteration_number = round(iteration_number / Qc);
% only report average iteration number over all sub-problems
% average per iteration per sub-cube
aruntime = aruntime/iteration_number;
atime_facet = atime_facet/iteration_number;
atime_data = atime_data/iteration_number;
%
vruntime = (sum_runtime_sqr - iteration_number*aruntime^2)/(iteration_number - 1);
vtime_facet = (sum_facet_sqr - iteration_number*atime_facet^2)/(iteration_number - 1);
vtime_data = (sum_data_sqr - iteration_number*atime_data^2)/(iteration_number - 1);
%
acpu_time = total_cpu_time/iteration_number;
vcpu_time = (sum_cpu_sqr - iteration_number*acpu_time^2)/(iteration_number - 1);

for k = 1:numel(Qc)
    fprintf("Qc = %i, iteration_number = %i \n", ...
        Qc, iteration_number/Qc)
    fprintf("aSTD=%1.2e, aKURT = %1.2e, aSKEW = %1.2e \n", ...
        astd_res, akurt_res, askew_res)
    fprintf(" aruntime (s) = %.2f, std_runtime (s) = %1.2e, acpu_time (s) = %.2f, std_cpu_time (s) = %1.2e \n", ...
        aruntime, sqrt(vruntime), acpu_time, sqrt(vcpu_time));
    fprintf(" total_runtime (h) = %2.2f, total_cpu_time (h) = %2.2f \n", ...
        total_runtime/3600, total_cpu_time/3600)
end

%% Saving results
save(metric_filename, '-v7.3', 'aruntime', 'vruntime', 'acpu_time', ...
    'vcpu_time', 'iteration_number', ...
    'atime_facet', 'atime_data', ...
    'vtime_facet', 'vtime_data', ...
    'total_runtime', 'total_cpu_time', 'std_res', 'kurt_res', 'skew_res');
end
