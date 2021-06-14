function [x, residual, asnr, ssnr, asnr_log, ssnr_log, acpu, scpu, arun, srun, total_cpu_time, total_runtime, iteration_number] = ...
    aggregate_results_spectral(algo, filename, ncores_data, ncores_prior, x0, squared_operator_norm, Qc, upsilon0)
%%
% Produce the images and metrics reported in the MNRAS paper
% ``A Faceted Prior for Scalable Wideband Imaging: Application to Radio
% Astronomy''
%
%-------------------------------------------------------------------------%
%%
% Author: P.-A. Thouvenin.
% Last modified: [../../2020]
%-------------------------------------------------------------------------%
%%
% NOTES:
%
%
% TODO: load operator norm, compute snr (print in a .txt file), compute
% timing, relate to number of CPUs
%
%-------------------------------------------------------------------------%
%%
N = [size(x0, 1), size(x0, 2)];
nChannels = size(x0, 3);
x = zeros(N(1), N(2), nChannels);
residual = zeros(N(1), N(2), nChannels);

% Reconstruction metrics
% DR = @(x, res, n) sqrt(prod(N))*squeeze(max(max(x,[],2),[],1))*n./norm2D(res); % per band input
norm2D = @(x) squeeze(sqrt(sum(sum(x.^2, 2), 1)));
SNR = @(x, x0) 20*log10(norm2D(x0)./norm2D(x - x0));
% SNR_log = @(x, x0) 20*log10(norm2D(log10(x0+eps))./norm2D(log10(x+eps)-log10(x0+eps)));
SNR_log = @(x, x0) 20*log10(norm2D(log10(1 + x0./reshape(upsilon0, [1, 1, numel(upsilon0)])))./norm2D(log10(1 + x./reshape(upsilon0, [1, 1, numel(upsilon0)]))-log10(1 + x0./reshape(upsilon0, [1, 1, numel(upsilon0)]))));

% asnr = zeros(numel(Qx), 1);
% asnr_log = zeros(numel(Qx), 1);
% vsnr = zeros(numel(Qx), 1);
% vsnr_log = zeros(numel(Qx), 1);
% iteration_number = zeros(numel(Qx), 1);
% total_runtime = zeros(numel(Qx), 1);  % total runtime (in s)
% total_cpu_time = zeros(numel(Qx), 1); % total CPU time (in s)
% aruntime = zeros(numel(Qx), 1);       % average runtime (per iteration, in s)
% acpu_time = zeros(numel(Qx), 1);      % average cpu time (per iter., in s)
% vcpu_time = zeros(numel(Qx), 1);      % variance cpu time
% atime_facet = zeros(numel(Qx), 1);    % average time (facet, per iter., in s)
% atime_data = zeros(numel(Qx), 1);     % variance
% vtime_facet = zeros(numel(Qx), 1);    % average time (data, per iter., in s)
% vtime_data = zeros(numel(Qx), 1);     % variance

switch algo
    case 'sara'
    
        total_runtime = 0;  % assuming all sub-cubes are running in parallel
        total_cpu_time = 0; %
        arun = 0;

        atime_l11 = 0; % average time (per iteration)
        atime_data = 0;
        atime_master = 0;
        iteration_number = 0;

        sum_runtime_sqr = 0;
        sum_l11_sqr = 0;
        sum_master_sqr = 0;
        sum_data_sqr = 0;
        sum_cpu_sqr = 0;

        for l = 1:nChannels
            load(filename(l), 'xsol', 't_l11', 't_master', 't_data', 'end_iter', 'res')
            x(:,:,l) = xsol;
            residual(:,:,l) = res/squared_operator_norm(l);
            
            % compute mean and variance over all the files? just 1 for now
            arun = arun + sum(end_iter(end_iter > 0)); % average runtime per iteration, over all sub-problems
            sum_runtime_sqr = sum_runtime_sqr + sum(end_iter(end_iter > 0).^2);
            
            atime_l11 = atime_l11 + sum(t_l11(t_l11 > 0));
            atime_master = atime_master + sum(t_master(t_master > 0));
            atime_data = atime_data + sum(t_data(t_data > 0));
            sum_l11_sqr = sum_l11_sqr + (ncores_prior^2)*sum(t_l11(t_l11 > 0).^2);
            sum_master_sqr = sum_master_sqr + sum(t_master(t_master > 0).^2);
            sum_data_sqr = sum_data_sqr + sum(t_data(t_data > 0).^2);
            
            % average number of iterations over all sub-problems
            iteration_number = iteration_number + sum(end_iter > 0);
            total_runtime = max(total_runtime, sum(end_iter(end_iter > 0)));
            total_cpu_time = total_cpu_time + ncores_prior*sum(t_l11(t_l11 > 0)) + ...
                + sum(t_master(t_master > 0)) ...
                + sum(t_data(t_data > 0));
            sum_cpu_sqr = sum_cpu_sqr + sum((ncores_prior*t_l11(t_l11 > 0) + t_master(t_master > 0) ...
                + t_data(t_data > 0)).^2);
        end
        
        arun = arun/iteration_number;
        % atime_l11 = atime_l11/iteration_number;
        % atime_master = atime_master/iteration_number;
        % atime_data = atime_data/iteration_number;
        %
        srun = sqrt((sum_runtime_sqr - iteration_number*arun^2)/(iteration_number - 1));
        % stime_l11 = sqrt((sum_l11_sqr - iteration_number*atime_l11^2)/(iteration_number - 1));
        % stime_master = sqrt((sum_master_sqr - iteration_number*atime_master^2)/(iteration_number - 1));
        % stime_data = sqrt((sum_data_sqr - iteration_number*atime_data^2)/(iteration_number - 1));
        %
        acpu = total_cpu_time/iteration_number;
        scpu = sqrt((sum_cpu_sqr - iteration_number*acpu^2)/(iteration_number - 1));
        iteration_number = round(iteration_number/nChannels); % average number of iterations for each SARA problem
        
        % compute SNR
        a = SNR(x, x0);
        asnr = mean(a);
        ssnr = std(a);
        a = SNR_log(x, x0);
        asnr_log = mean(a);
        ssnr_log = std(a);

        ncpu = nChannels*(ncores_data + ncores_prior + 1);

    case 'hypersara'

        total_runtime = 0;  % assuming all sub-cubes are running in parallel
        total_cpu_time = 0; %
        arun = 0;

        atime_l21 = 0; % average time (per iteration)
        atime_nuclear = 0;
        atime_data = 0;
        atime_master = 0;
        iteration_number = 0;

        sum_runtime_sqr = 0;
        sum_l21_sqr = 0;
        sum_nuclear_sqr = 0;
        sum_master_sqr = 0;
        sum_data_sqr = 0;
        sum_cpu_sqr = 0;

        id = split_range_interleaved(Qc(c), nChannels);
        x = zeros(size(x0));
        residual = zeros(size(x0));

        for ind = 1:numel(id)
            load(filename(ind), 'xsol', 'res', 't_l21', 't_nuclear', 't_master', 't_data', 'end_iter')
            x(:,:,id{ind}) = xsol;
            residual(:,:,id{ind}) = res;
            clear xsol res;
        
            % compute mean and variance over all the files
            arun = arun + sum(end_iter(end_iter > 0)); % average runtime per iteration, over all sub-problems
            sum_runtime_sqr = sum_runtime_sqr + sum(end_iter(end_iter > 0).^2);
            
            atime_l21 = atime_l21 + sum(t_l21(t_l21 > 0));
            atime_nuclear = atime_nuclear + sum(t_nuclear(t_nuclear > 0));
            atime_master = atime_master + sum(t_master(t_master > 0));
            atime_data = atime_data + sum(t_data(t_data > 0));
            sum_l21_sqr = sum_l21_sqr + sum(t_l21(t_l21 > 0).^2);
            sum_nuclear_sqr = sum_nuclear_sqr + sum(t_nuclear(t_nuclear > 0).^2);
            sum_master_sqr = sum_master_sqr + sum(t_master(t_master > 0).^2);
            sum_data_sqr = sum_data_sqr + (ncores_data^2)*sum(t_data(t_data > 0).^2); % results from data cores have already been averaged over data cores -> need to be taken into account in the variances
            
            % average number of iterations over all sub-problems
            iteration_number = iteration_number + sum(end_iter > 0);
            total_runtime = max(total_runtime, sum(end_iter(end_iter > 0)));
            total_cpu_time = total_cpu_time + sum(t_l21(t_l21 > 0)) + ...
                sum(t_nuclear(t_nuclear > 0)) + sum(t_master(t_master > 0)) ...
                + ncores_data*sum(t_data(t_data > 0));
            sum_cpu_sqr = sum_cpu_sqr + sum((t_l21(t_l21 > 0) + t_nuclear(t_nuclear > 0) + t_master(t_master > 0) ...
                + ncores_data*t_data(t_data > 0)).^2);

        end

        % iteration_number(k) = round(iteration_number(k) / Qc(k)); % only report average iteration number over all sub-problems
        arun = arun/iteration_number;
        % atime_l21 = atime_l21/iteration_number;
        % atime_nuclear = atime_nuclear/iteration_number;
        % atime_master = atime_master/iteration_number;
        % atime_data = atime_data/iteration_number;
        %
        srun = sqrt((sum_runtime_sqr - iteration_number*arun^2)/(iteration_number - 1));
        % stime_l21 = sqrt((sum_l21_sqr - iteration_number*atime_l21^2)/(iteration_number - 1));
        % stime_nuclear = sqrt((sum_nuclear_sqr - iteration_number*atime_nuclear^2)/(iteration_number - 1));
        % stime_master = sqrt((sum_master_sqr - iteration_number*atime_master^2)/(iteration_number - 1));
        % stime_data = sqrt((sum_data_sqr - iteration_number*atime_data^2)/(iteration_number - 1));  
        % 
        acpu = total_cpu_time/iteration_number;
        scpu = sqrt((sum_cpu_sqr - iteration_number*acpu^2)/(iteration_number - 1));
        iteration_number = round(iteration_number/nChannels); % average number of iterations for each problem

        ncpu = min(ncores_prior + ncores_data + 1, 36);

        % compute SNR
        a = SNR(x, x0);
        asnr = mean(a);
        ssnr = std(a);
        a = SNR_log(x, x0);
        asnr_log = mean(a);
        ssnr_log = std(a);

        residual= res./reshape(squared_operator_norm, [1, 1, nChannels]);
end

%% Display results (table)
fprintf("%s, asnr = %2.2f (%1.2e), asnr_log = %2.2f (%1.2e), cpu_cores = %i, iteration_number = %i \n", ...
    algo, asnr, ssnr, asnr_log, ssnr_log, ncpu, iteration_number)
fprintf("arun (s) = %.2f (%1.2e), total_runtime (h) = %2.2f \n", ...   
    arun, srun, total_runtime/3600);
fprintf("acpu (s) = %.2f (%1.2e) , total_cpu_time (h) = %2.2f \n", ...
    acpu, scpu, total_cpu_time/3600);

end